function [optim_results, grid_points, grid_scores, gauss_fit, metadata] = fitModelQRG(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELQRG Strategy for fitting inference models by inverse-CDF sampling from the prior
%using a quasi-random grid (QRG), then maximizing a gaussian-process fit to those points.

fields = fieldnames(distribs);

% Get a UID for the model being fit, signals, and choices so that we can restart from checkpoints
if iscell(signals)
    allsigs = vertcat(signals{:});
    allchoices = vertcat(choices{:});
    input_id = string2hash([base_params(1).model, strjoin(fields), num2str([allsigs(:)' allchoices'])]);
    nTrials = length(allchoices);
else
    input_id = string2hash([base_params(1).model, strjoin(fields), num2str([signals(:)' choices'])]);
    nTrials = length(choices);
end
chkpt = fullfile('sample-checkpoints', sprintf('%x-qrg.mat', input_id));

%% Set up evaluation, handling stochastic vs deterministic case
if Model.isStochastic(base_params)
    disp('fitModelQRG (stochastic case)');
    ibs_repeats = 5;
    
    eval_fun = @Fitting.choiceModelLogProbIBS;
    eval_args = {[], ibs_repeats};
    final_eval_args = {[], round(10*sqrt(nTrials))};
else
    disp('fitModelQRG (deterministic case)');
    eval_fun = @Fitting.choiceModelLogProb;
    eval_args = {1};
    final_eval_args = {1};
end

%% Evaluate parameters on a quasi-random grid (QRG)

% To gracefully save/load from partially computed results, we initialize to NaN values, (maybe) load
% from a file, then only evaluate those rows that are still NaN, with periodic checkpointing.
if exist(chkpt, 'file')
    ld = load(chkpt);
    grid_points = ld.grid_points;
    logpdf = ld.logpdf;
    loglike = ld.loglike;
    variance = ld.variance;
    
    fprintf('fitModelQRG :: %s\tloaded %d of %d evaluations\n', chkpt, sum(~isnan(loglike)), length(loglike));
else
    % Start with a QRG on the unit hypercube
    hyperQRG = lhsdesign(10000, length(fields));
    
    % Use points inside hypercube to get grid on actual parameters via inverse CDF of each field's prior
    grid_points = zeros(size(hyperQRG));
    for iF=1:length(fields)
        grid_points(:,iF) = distribs.(fields{iF}).priorinv(hyperQRG(:,iF));
    end
    
    % Initialize all evaluations to NaN
    logpdf = nan(size(hyperQRG,1), 1);
    loglike = nan(size(hyperQRG,1), 1);
    variance = nan(size(hyperQRG,1), 1);
end

tstart = tic;
needs_eval = find(isnan(loglike));
for ii=1:length(needs_eval)
    iVal = needs_eval(ii);
    [logpdf(iVal), loglike(iVal), variance(iVal)] = ...
        eval_fun(Fitting.setParamsFields(base_params, fields, grid_points(iVal, :)), distribs, signals, choices, eval_args{:});
    
    % Periodically print diagnostics and save things
    if mod(ii, 100) == 0 || ii == length(needs_eval)
        telapse = toc(tstart);
        fprintf('%s :: eval %d of %d\tETA=%.1fs\n', chkpt, ii, length(needs_eval), (length(needs_eval)-ii)*telapse/ii);
        save(chkpt, 'grid_points', 'logpdf', 'loglike', 'variance');
    end
end

% Collect results into 'grid_scores' struct
grid_scores.logpdf = logpdf;
grid_scores.loglike = loglike;
grid_scores.variance = variance;
clear logpdf loglike variance;

%% Find best grid point based on existing (rough) scores
disp('fitModelQRG: search grid');

[~, sort_ll] = sort(grid_scores.loglike, 'descend');
best_ll_idx = sort_ll(1);
optim_results.mle_params = Fitting.setParamsFields(base_params, fields, grid_points(best_ll_idx, :));

[~, sort_lp] = sort(grid_scores.logpdf, 'descend');
best_lp_idx = sort_lp(1);
optim_results.map_params = Fitting.setParamsFields(base_params, fields, grid_points(best_lp_idx, :));

%% Estimate moment-matched gaussian posterior
disp('fitModelQRG: Gauss fit');

relative_ll = grid_scores.loglike - max(grid_scores.loglike);
% Equation (41) in "Unbiased and Efficient Log-Likelihood Estimation with Inverse Binomial Sampling"
est_likelihood = exp(relative_ll + 1/2*grid_scores.variance);
weight = est_likelihood / sum(est_likelihood);

% Since 'grid_points' are from the prior, we construct an estimator of posterior moments by
% weighting each point by an estimate of the likelihood
mu = weight' * grid_points;

% Subtract mean and compute 2nd moment (covariance)
grid0 = grid_points - mu;

gauss_fit.mu = mu;
gauss_fit.Sigma = squeeze(sum(weight .* grid0 .* reshape(grid0, [], 1, length(fields))));

optim_results.mean_params = Fitting.setParamsFields(base_params, fields, mu);

%% Prepare GP stuff
disp('fitModelQRG: GP Prep');

% Set initial length scales as 1/4th of the range of plausible bounds for each parameter (this was
% chosen after fittting the GP many times and observing that 1/4 was a good estimate of where it
% would converge to anyway).
init_scale = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/4, fields)';
% Set convergence tolerance as 1/1000th of the plausible range
tolerance = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/1000, fields)';

observation_variance = max(.01, dot(weight, grid_scores.variance));

% The GP fit shouldn't waste time trying to get the tails right, and the defaults for fitrgp only
% select 2000 points anyway. So, only bother trying to fit the top 2000 scoring points.
gp_idx = sort_ll(1:2000);

% The basis (mean) function for the GP will be a quadratic model, H*beta, where H is [1 x x.^2].
% (This is what happens inside fitrgp with the 'PureQuadratic' basis flag). Here we guess beta
% according to the gaussian parameters estimate above, then refine it by least squares fitting
H = [ones(length(gp_idx),1) grid_points(gp_idx, :) grid_points(gp_idx, :).^2];
% The following reparameterizes -1/2*(x-mu)^2/sigma^2 into H*beta format.
gauss_beta = [max(grid_scores.loglike) - sum(gauss_fit.mu'.^2./diag(gauss_fit.Sigma))/2;
    gauss_fit.mu'./diag(gauss_fit.Sigma);
    -1./(2*diag(gauss_fit.Sigma))];
% Optimize 'beta' parameter with constraint that quadratic terms are negative (enforcing bounded
% (gaussian) probability.
lb = -inf(size(gauss_beta));
ub = inf(size(gauss_beta));
ub(length(fields)+2:end) = 0;
opts = optimoptions('fmincon', 'display', 'none');
fit_beta = fmincon(@(beta) sum((grid_scores.loglike(gp_idx) - H*beta).^2), gauss_beta, ...
    [], [], [], [], lb, ub, [], opts);

% Estimate how much variance to ascribe to each term in the GP. Majority of the variance of 'logpdf'
% and 'loglike' themselves comes from the quadratic basis, which we first subtract out to get the
% 'residual' variance that the GP ought to explain.
residual_variance = max(.01, var(grid_scores.loglike(gp_idx) - H*fit_beta) - observation_variance);

% Divide up the unexplained variance per each parameter's marginal.
init_sigma = sqrt(residual_variance / length(fields)) * ones(1, length(fields));

% Create initial kernel params: vertical stack alternating scale/sigma once per field, then at the
% final index 2*|fields|+1, the 'sigma' corresponding to the joint prob fluctuations that aren't
% explained by the sum of marginals.
init_gp_kernel_params = [init_scale(:)'; init_sigma];
init_gp_kernel_params = [init_gp_kernel_params(:); init_sigma(1)];
% Kernel expects *log* scales and *log* sigmas
init_gp_kernel_params = log(init_gp_kernel_params);

gp_args = {'Sigma', sqrt(observation_variance), 'ConstantSigma', true, 'PredictMethod', 'exact', ...
    'KernelFunction', @Fitting.marginalJointKernel, 'KernelParameters', init_gp_kernel_params, ...
    'Basis', 'PureQuadratic', 'Beta', fit_beta, 'FitMethod', 'none'};

%% Optimize GP kernel

% The GP has 3 distinct sets of hyperparameters: scalar 'Sigma' or observation noise, vector 'Beta'
% parameterizing the mean or basis function, and the kernel hyperparameters. At this point, beta has
% been fit by least squares, and 'Sigma' is set according to grid_scores.variance.  The GP scale is
% set according to the prior bounds and will not be messed-with. These should not be changed. The
% final parameters to optimize are the variances of the GP series themselves, i.e. the kernel
% sigmas. Since built-in fitrgp either optimizes all hyperparameters or none of them, we do our own
% search here ot optimize only the sigmas.

% Indices into init_gp_kernel_params containing log(sigma)s (see @Fitting.marginalJointKernel)
sigidx = [2:2:length(init_gp_kernel_params) length(init_gp_kernel_params)];

% Initially search a range from 0.25x to 4x the initial guess. This will shrink each iteration.
scale = 2.^(-2:.25:2);

% Get index into the big 'gp_args' cell array containing init_gp_kernel_params
kernel_arg_idx = find(strcmpi(gp_args, 'KernelParameters'))+1;

% Keep a record of how parameters changed each iteration (debugging)
record(1, :) = init_gp_kernel_params;

for iter=1:20
    loss = zeros(length(sigidx), length(scale));
    for isig=1:length(sigidx)
        for iscale=1:length(scale)
            % Set the value of sigma(isig) scaled by scale(iscale) (parameter is log(sigma))
            this_para = init_gp_kernel_params;
            this_para(sigidx(isig)) = this_para(sigidx(isig)) + log(scale(iscale));
            
            % Construct a new random train/test split (90%/10%)
            idxtrain = randperm(length(gp_idx), round(0.9*length(gp_idx)));
            idxtest = setdiff(1:length(gp_idx), idxtrain);
            
            % Create ('train') this GP
            gp_args{kernel_arg_idx} = this_para;
            gp_ll_xv = fitrgp(grid_points(gp_idx(idxtrain), :), grid_scores.loglike(gp_idx(idxtrain)), gp_args{:});
            
            % Compute squared error on held out points
            loss(isig,iscale) = sum((gp_ll_xv.predict(grid_points(gp_idx(idxtest), :)) - grid_scores.loglike(gp_idx(idxtest))).^2);
        end
        % Apply some smoothing to the loss to avoid chasing spurious minima
        loss(isig, :) = smooth(loss(isig, :), 5);
    end
    % Select best 'scale' along each dimension
    [~,ibest] = min(loss, [], 2);
    init_gp_kernel_params(sigidx) = init_gp_kernel_params(sigidx) + log(scale(ibest))';
    scale = exp(log(scale) * .9);
    
    record(iter+1, :) = init_gp_kernel_params;
    clf;
    subplot(2,1,1);
    plot(record(:,sigidx), 'marker', '.');
    xlabel('iteration');
    ylabel('log(\sigma)');
    legend([fields' {'joint'}], 'location', 'best');
    subplot(2,1,2);
    imagesc(loss); axis image; colorbar;
    drawnow;
end

keyboard;

% Set args to use the best sigmas
gp_args{kernel_arg_idx} = init_gp_kernel_params;

%% Search LL landscape with GP
disp('fitModelQRG: GP Fit [LL]');

% Fit Gaussian process to all log likelihood evaluations
gp_ll = fitrgp(grid_points(gp_idx, :), grid_scores.loglike(gp_idx), gp_args{:});

% Each iteration, restart the search from the 10 best samples from above
searchpoints = grid_points(sort_ll(1:10), :);
nsearchjitter = 10;

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
lb = cellfun(@(f) distribs.(f).lb, fields);
ub = cellfun(@(f) distribs.(f).ub, fields);
gp_mle(1,:) = gpmax(gp_ll, searchpoints, nsearchjitter, [], [], [], [], lb, ub, [], opts);

disp('fitModelQRG: GP Search [LL]');

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_ll.predict(gp_mle(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_mle=%.1f+/-%.1f', itr, max(abs(delta)), ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [~, new_ll, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_mle(itr,:)), ...
        distribs, signals, choices, eval_args{:});

    gp_ll = updateGPModel(gp_ll, gp_mle(itr,:), new_ll);
    
    fprintf('\tactual_mle=%.1f+/-%.1f\n', new_ll, sqrt(new_var));
    
    % Search again for the new best point
    gp_mle(itr+1,:) = gpmax(gp_ll, searchpoints, nsearchjitter, [], [], [], [], lb, ub, [], opts);
    
    delta = gp_mle(itr+1,:) - gp_mle(itr,:);
    itr = itr+1;
end

optim_results.gp_mle_params = Fitting.setParamsFields(base_params, fields, gp_mle(end,:));

metadata.gp_ll = gp_ll;

%% Repeat the above for the MAP
disp('fitModelQRG: GP Fit [MAP]');

% Fit Gaussian process to all log likelihood evaluations
gp_lp = fitrgp(grid_points(gp_idx,:), grid_scores.logpdf(gp_idx), gp_args{:});

% Each iteration, restart the search from the 10 best samples from above
searchpoints = grid_points(sort_lp(1:10), :);
nsearchjitter = 10;

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
gp_map(1,:) = gpmax(gp_ll, searchpoints, nsearchjitter, [], [], [], [], lb, ub, [], opts);

disp('fitModelQRG: GP Search [MAP]');

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_lp.predict(gp_map(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_logp=%.1f+/-%.1f', itr, max(abs(delta)), ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [new_lp, ~, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_map(itr,:)), ...
        distribs, signals, choices, eval_args{:});

    gp_lp = updateGPModel(gp_lp, gp_map(itr,:), new_lp);
    
    fprintf('\tactual_logp=%.1f+/-%.1f\n', new_lp, sqrt(new_var));
    
    % Search again for the new best point
    gp_map(itr+1,:) = gpmax(gp_ll, searchpoints, nsearchjitter, [], [], [], [], lb, ub, [], opts);
    
    delta = gp_map(itr+1,:) - gp_map(itr,:);
    itr = itr+1;
end

optim_results.gp_map_params = Fitting.setParamsFields(base_params, fields, gp_map(end,:));

metadata.gp_map = gp_map;

%% Re-run and store final evaluations for each model
fit_names = fieldnames(optim_results);
for iFit=1:length(fit_names)
    fprintf('fitModelQRG: final evals [%s]\n', fit_names{iFit});
    [lp, ll, ll_var] = eval_fun(optim_results.(fit_names{iFit}), distribs, signals, choices, final_eval_args{:});
    for iP=1:length(optim_results.(fit_names{iFit}))
        optim_results.(fit_names{iFit})(iP).lp = lp;
        optim_results.(fit_names{iFit})(iP).ll = ll;
        optim_results.(fit_names{iFit})(iP).ll_var = ll_var;
    end
end

end

function [bestx, bestval] = gpmax(gp, startpts, jitters, varargin)
% Helper to search around a set of seed points to find the maximum of a GP
scale = exp(gp.KernelInformation.KernelParameters(1:2:end-2))';

xval = zeros(size(startpts,1)*(1+jitters), size(startpts, 2));
fval = zeros(size(xval,1),1);

idx = 1;
for istart=1:size(startpts,1)
    % Start a search exactly at the ith 'seed' point
    [xval(idx,:), fval(idx,:)] = fmincon(@(x) -gp.predict(x), startpts(istart,:), varargin{:});
    idx = idx+1;

    % Jitter around this point 'jitters' number of times and start a search at each one
    for ijit=1:jitters
        [xval(idx,:), fval(idx,:)] = fmincon(@(x) -gp.predict(x), startpts(istart,:)+randn(size(scale)).*scale, varargin{:});
        idx = idx+1;
    end
end
[bestval, idxbest] = min(fval);
bestx = xval(idxbest,:);
bestval = -bestval;
end

function mdl = updateGPModel(mdl,Xadd,Yadd)
%UPDATEGPMoDeL Helper to 'add' a data point to a GP model without updating hyperparameters. Thanks to
%https://www.mathworks.com/matlabcentral/answers/424553-how-to-add-points-to-a-trained-gp
kernelfun = mdl.KernelFunction;
kernelparams = mdl.KernelInformation.KernelParameters;
sigma = mdl.Sigma;
basis = mdl.BasisFunction;
beta = mdl.Beta;
Xall = [mdl.X; Xadd];
Yall = [mdl.Y; Yadd];

mdl = fitrgp(Xall, Yall, 'Sigma', sigma, 'Basis', basis, 'Beta', beta, 'KernelParameters', kernelparams, ...
    'KernelFunction', kernelfun, 'ConstantSigma', true, 'PredictMethod', 'exact', 'FitMethod', 'none');
end