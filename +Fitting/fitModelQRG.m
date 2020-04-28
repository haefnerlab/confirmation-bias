function [optim_results, grid_points, grid_scores, metadata] = fitModelQRG(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELQRG Strategy for fitting inference models by inverse-CDF sampling from the prior
%using a quasi-random grid (QRG), then maximizing a gaussian-process fit to those points.

fields = fieldnames(distribs);
nF = length(fields);

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
    hyperQRG = lhsdesign(10000, nF);
    
    % Use points inside hypercube to get grid on actual parameters via inverse CDF of each field's prior
    grid_points = zeros(size(hyperQRG));
    for iF=1:nF
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
disp('fitModelQRG: Quad fit');

% Fits shouldn't waste time trying to get the tails right, and the defaults for fitrgp only select
% 2000 points anyway. So, only bother trying to fit the top 2000 scoring points.
fit_idx = sort_ll(1:2000);

% Fit a gaussian, i.e. fit a quadratic function to the log-likelihood. ll≈H*beta, where H is a
% quadratic basis.
x0 = ones(length(fit_idx), 1);
x1 = grid_points(fit_idx,:);
x2 = grid_points(fit_idx,:).^2;
% Construct regressors H
H = [x0 x1 x2];

% Set up constrained least-squares fit to the log likelihood, constrained so that the curvature is
% downward (ensuring probability under the implied gaussian is bounded).
lb = -inf(size(H, 2), 1);
ub = inf(size(lb));
ub(nF+2:end) = 0; % the coefficients on the variances (a*x_i^2) must be negative (a<0)
init_beta = [mean(grid_scores.loglike(fit_idx));
    zeros(2*nF,1)];
init_beta(ub==0) = -1;

% Use 'fit_weight' to encourage fit to focus on the higher LL points, but only slightly.
fit_weight = 1./(1+exp(-zscore(grid_scores.logpdf(fit_idx))));
fit_weight = fit_weight / sum(fit_weight);

% Do optimization
opts = optimoptions('fmincon', 'display', 'final', 'MaxFunctionEvaluations', inf);
fit_beta = fmincon(@(beta) sum(fit_weight.*(grid_scores.loglike(fit_idx) - H*beta).^2), init_beta, ...
    [], [], [], [], lb, ub, [], opts);

% Create local function x -> baseline replicating H*beta above
    function b = quadBasis(x)
        b = [ones(size(x,1),1) x x.^2] * fit_beta;
    end

% Use best-fit quadratic function as the basis of the GP fit... implemented by subtracting off the
% quadratic prediction then fitting a constant-baseline GP to the residuals
residual_loglike = grid_scores.loglike(fit_idx) - quadBasis(grid_points(fit_idx, :));

%% Prepare GP stuff
disp('fitModelQRG: GP Prep');

% Set initial length scales as 1/4th of the range of plausible bounds for each parameter (this was
% chosen after fittting the GP many times and observing that 1/4 was a good estimate of where it
% would converge to anyway).
init_scale = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/4, fields)';
% Set convergence tolerance as 1/1000th of the plausible range
tolerance = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/1000, fields)';

observation_variance = max(1, dot(fit_weight, grid_scores.variance(fit_idx)));

% Estimate how much variance to ascribe to each term in the GP. Majority of the variance of 'logpdf'
% and 'loglike' themselves comes from the quadratic basis which was subtracted out above.
residual_variance = max(1, var(residual_loglike) - observation_variance);

% Divide up the unexplained variance per each parameter's marginal.
init_sigma = sqrt(residual_variance / nF) * ones(1, nF);

% Create initial kernel params: vertical stack alternating scale/sigma once per field, then at the
% final index 2*|fields|+1, the 'sigma' corresponding to the joint prob fluctuations that aren't
% explained by the sum of marginals.
init_gp_kernel_params = [init_scale(:)'; init_sigma];
init_gp_kernel_params = [init_gp_kernel_params(:); init_sigma(1)];
% Kernel expects *log* scales and *log* sigmas
init_gp_kernel_params = log(init_gp_kernel_params);

gp_args = {'Sigma', sqrt(observation_variance), 'ConstantSigma', true, 'PredictMethod', 'exact', ...
    'KernelFunction', @Fitting.marginalJointKernel, 'KernelParameters', init_gp_kernel_params, ...
    'Basis', 'constant', 'FitMethod', 'exact', 'CrossVal', 'on'};

%% Search LL landscape with GP
disp('fitModelQRG: GP Fit');

% Fit Gaussian process to all log likelihood evaluations, cross validated (default is 10-fold)
xv_gp_ll = fitrgp(grid_points(fit_idx, :), residual_loglike, gp_args{:});

% Aggregate kernel hyperparameters: take their geometric mean
avg_gp_kernel_params = zeros(size(init_gp_kernel_params));
nFold = length(xv_gp_ll.Trained);
for iFold=1:nFold
    avg_gp_kernel_params = avg_gp_kernel_params + xv_gp_ll.Trained{iFold}.KernelInformation.KernelParameters / nFold;
end

% Construct actual GP object that will be used for searching, with further hyperparameter tuning
% turned off
gp_args{find(strcmp(gp_args, 'KernelParameters'))+1} = avg_gp_kernel_params;
gp_args{find(strcmp(gp_args, 'FitMethod'))+1} = 'none';
gp_args{find(strcmp(gp_args, 'CrossVal'))+1} = 'off';
gp_ll = fitrgp(grid_points(fit_idx, :), residual_loglike, gp_args{:});

% Each iteration, restart the search from the 10 best samples from above
searchpoints = grid_points(sort_ll(1:10), :);
nsearchjitter = 10;

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
jitterscale = exp(gp_ll.KernelInformation.KernelParameters(1:2:end-2))';
lb = cellfun(@(f) distribs.(f).lb, fields);
ub = cellfun(@(f) distribs.(f).ub, fields);
gp_mle(1,:) = fminconrepeats(@(x) -(quadBasis(x)+gp_ll.predict(x)), searchpoints, nsearchjitter, jitterscale, ...
    [], [], [], [], lb, ub, [], opts);

disp('fitModelQRG: GP Search [LL]');

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_ll.predict(gp_mle(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_ll=%.1f+/-%.1f', itr, max(abs(delta)), quadBasis(gp_mle(itr,:))+ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [~, new_ll, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_mle(itr,:)), ...
        distribs, signals, choices, eval_args{:});

    gp_ll = updateGPModel(gp_ll, gp_mle(itr,:), new_ll-quadBasis(gp_mle(itr, :)));
    
    dists = sum(((grid_points-gp_mle(itr,:))./exp(gp_ll.KernelInformation.KernelParameters(1:2:end-2))').^2,2);
    [~,iclosest] = min(dists);
    
    fprintf('\tactual_ll=%.1f+/-%.1f\tclosest_ll=%.1f+/-%.1f\n', new_ll, sqrt(new_var), grid_scores.loglike(iclosest), sqrt(grid_scores.variance(iclosest)));
    
    disp([gp_mle(itr,:); grid_points(iclosest,:)]);
    
    % Search again for the new best point
    gp_mle(itr+1,:) = fminconrepeats(@(x) -(quadBasis(x)+gp_ll.predict(x)), searchpoints, nsearchjitter, jitterscale, ...
        [], [], [], [], lb, ub, [], opts);    
    delta = gp_mle(itr+1,:) - gp_mle(itr,:);
    itr = itr+1;
end

optim_results.gp_mle_params = Fitting.setParamsFields(base_params, fields, gp_mle(end,:));

metadata.gp_ll = gp_ll;

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

function [bestx, bestval] = fminconrepeats(fcn, startpts, jitters, jitterscale, varargin)
% Helper to search around a set of seed points to find the maximum of a function
xval = zeros(size(startpts,1)*(1+jitters), size(startpts, 2));
fval = zeros(size(xval,1),1);

idx = 1;
for istart=1:size(startpts,1)
    % Start a search exactly at the ith 'seed' point
    [xval(idx,:), fval(idx,:)] = fmincon(fcn, startpts(istart,:), varargin{:});
    idx = idx+1;

    % Jitter around this point 'jitters' number of times and start a search at each one
    for ijit=1:jitters
        [xval(idx,:), fval(idx,:)] = fmincon(fcn, startpts(istart,:)+randn(size(jitterscale)).*jitterscale, varargin{:});
        idx = idx+1;
    end
end
[bestval, idxbest] = min(fval);
bestx = xval(idxbest,:);
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