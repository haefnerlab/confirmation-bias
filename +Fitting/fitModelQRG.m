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

% Set initial length scales as 1/20th of the range of plausible bounds for each parameter
init_scale = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/20, fields)';
% Set convergence tolerance as 1/1000th of the plausible range
tolerance = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/1000, fields)';

% The GP fit shouldn't waste time trying to get the tails right, and the defaults for fitrgp only
% select 2000 points anyway. So, only bother trying to fit the top 2000 scoring points (plus keep
% around 100 of them for cross validation).
gp_idx = sort_ll(1:2100);
gp_holdout = 100/2100;

% Estimate how much variance to ascribe to each term in the GP. Majority of the variance of 'logpdf'
% and 'loglike' themselves comes from the quadratic basis, which we can estimate by directly fitting
% a quadratic model to the log probabilties.
mdl = fitlm(grid_points(gp_idx,:), grid_scores.loglike(gp_idx), 'purequadratic');
residual_variance = var(grid_scores.loglike(gp_idx) - mdl.predict());

% Divide up the unexplained variance per each parameter's marginal
init_sigma = sqrt(residual_variance / length(fields)) * ones(1, length(fields));

% Create initial kernel params: vertical stack alternating scale/sigma once per field, then at the
% final index 2*|fields|+1, the 'sigma' corresponding to the joint prob fluctuations that aren't
% explained by the sum of marginals.
init_gp_kernel_params = [init_scale(:)'; init_sigma];
init_gp_kernel_params = [init_gp_kernel_params(:); init_sigma(1)];
% Kernel expects *log* scales and *log* sigmas
init_gp_kernel_params = log(init_gp_kernel_params);

gp_args = {'Sigma', max(0.01, sqrt(dot(weight, grid_scores.variance))), 'ConstantSigma', true, 'PredictMethod', 'exact', ...
    'Basis', 'PureQuadratic', 'Beta', mdl.Coefficients.Estimate, 'HoldOut', gp_holdout, 'KernelFunction', @Fitting.marginalJointKernel, ...
    'KernelParameters', init_gp_kernel_params};

%% Search LL landscape with GP
disp('fitModelQRG: GP Fit [LL]');

% Fit Gaussian process to all log likelihood evaluations
gp_ll = fitrgp(grid_points(gp_idx, :), grid_scores.loglike(gp_idx), gp_args{:});
gp_ll = gp_ll.Trained{1};

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

    gp_ll = updateGPRMdl(gp_ll, gp_mle(itr,:), new_ll, mod(itr,20)==0);
    
    fprintf('\tactual_mle=%.1f+/-%.1f\n', new_ll, sqrt(new_var));
    
    % Search again for the new best point
    gp_mle(itr+1,:) = gpmax(gp_ll, searchpoints, nsearchjitter, [], [], [], [], lb, ub, [], opts);
    
    delta = gp_mle(itr+1,:) - gp_mle(itr,:);
    itr = itr+1;
end

optim_results.gp_mle_params = Fitting.setParamsFields(base_params, fields, gp_mle(end,:));

%% Repeat the above for the MAP
disp('fitModelQRG: GP Fit [MAP]');

% Fit Gaussian process to all log likelihood evaluations
gp_lp = fitrgp(grid_points(gp_idx,:), grid_scores.logpdf(gp_idx), gp_args{:});
gp_lp = gp_lp.Trained{1};

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

    gp_lp = updateGPRMdl(gp_lp, gp_map(itr,:), new_lp, mod(itr,20)==0);
    
    fprintf('\tactual_logp=%.1f+/-%.1f\n', new_lp, sqrt(new_var));
    
    % Search again for the new best point
    gp_map(itr+1,:) = gpmax(gp_ll, searchpoints, nsearchjitter, [], [], [], [], lb, ub, [], opts);
    
    delta = gp_map(itr+1,:) - gp_map(itr,:);
    itr = itr+1;
end

optim_results.gp_map_params = Fitting.setParamsFields(base_params, fields, gp_map(end,:));

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
scale = gp.KernelInformation.KernelParameters(1:end-1)';

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

function newmdl = updateGPRMdl(mdl,Xadd,Yadd,refitHypers)
%UPDATEGPRMDL Helper to 'add' a data point to a GP model without necessarily triggering a slow
%hyperparameter update. Thanks to
%https://www.mathworks.com/matlabcentral/answers/424553-how-to-add-points-to-a-trained-gp
if nargin < 4, refitHypers = true; end
kernelfun = mdl.KernelFunction;
kernelparams = mdl.KernelInformation.KernelParameters;
sigma = mdl.Sigma;
constsig = mdl.ModelParameters.ConstantSigma;
beta = mdl.Beta;
Xall = [mdl.X; Xadd];
Yall = [mdl.Y; Yadd];

if ~refitHypers
    fitarg = {'FitMethod', 'none'};
else
    % Let default fitting take care of it
    fitarg = {};
end

newmdl = fitrgp(Xall, Yall, 'Sigma', sigma, 'Beta', beta, 'KernelParameters', ...
    kernelparams, 'KernelFunction', kernelfun, 'ConstantSigma', constsig, ...
    'PredictMethod', 'exact', fitarg{:});
end