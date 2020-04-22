function [optim_results, grid_points, grid_scores, metadata] = fitModelQRG(base_params, signals, choices, distribs, metadata)
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

[~, best_ll_idx] = max(grid_scores.loglike);
optim_results.mle_params = Fitting.setParamsFields(base_params, fields, grid_points(best_ll_idx, :));

[~, best_lp_idx] = max(grid_scores.logpdf);
optim_results.map_params = Fitting.setParamsFields(base_params, fields, grid_points(best_lp_idx, :));

%% Search LL landscape with GP
disp('fitModelQRG: GP MLE');

% Set initial length scales as 1/20th of the range of plausible bounds for each parameter
init_scale = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/20, fields)';
% Set convergence tolerance as 1/1000th of the plausible range
tolerance = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/1000, fields)';

% Fit Gaussian process to all log likelihood evaluations
gp_ll = fitrgp(grid_points, grid_scores.loglike, 'KernelFunction', 'ardsquaredexponential', ...
    'KernelParameters', [init_scale std(grid_scores.loglike)], 'Sigma', max(0.01, mean(grid_scores.variance)));

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
lb = cellfun(@(f) distribs.(f).lb, fields);
ub = cellfun(@(f) distribs.(f).ub, fields);
gp_mle(1,:) = gpmax(gp_ll, grid_points(best_ll_idx, :), [], [], [], [], lb, ub, [], opts);

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_ll.predict(gp_mle(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_mle=%.1f+/-%.1f', itr, max(abs(delta)), ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [~, new_ll, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_mle(itr,:)), ...
        distribs, signals, choices, eval_args{:});
    refit = mod(itr,20)==0;
    gp_ll = updateGPRMdl(gp_ll, gp_mle(itr,:), new_ll, refit);
    
    fprintf('\tactual_mle=%.1f+/-%.1f\n', new_ll, sqrt(new_var));
    
    % Search again for the new best point
    gp_mle(itr+1,:) = gpmax(gp_ll, gp_mle(itr,:), [], [], [], [], lb, ub, [], opts);
    
    delta = gp_mle(itr+1,:) - gp_mle(itr,:);
    itr = itr+1;
end

optim_results.gp_mle_params = Fitting.setParamsFields(base_params, fields, gp_mle(end,:));

%% Repeat the above for the MAP
disp('fitModelQRG: GP MAP');

% Fit Gaussian process to all log likelihood evaluations
gp_lp = fitrgp(grid_points, grid_scores.logpdf, 'KernelFunction', 'ardsquaredexponential', ...
    'KernelParameters', [init_scale std(grid_scores.logpdf)], 'Sigma', max(0.01, mean(grid_scores.variance)));

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
gp_map(1,:) = gpmax(gp_lp, grid_points(best_lp_idx, :), [], [], [], [], lb, ub, [], opts);

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_lp.predict(gp_map(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_logp=%.1f+/-%.1f', itr, max(abs(delta)), ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [new_lp, ~, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_map(itr,:)), ...
        distribs, signals, choices, eval_args{:});
    refit = mod(itr,20)==0;
    gp_lp = updateGPRMdl(gp_lp, gp_map(itr,:), new_lp, refit);
    
    fprintf('\tactual_logp=%.1f+/-%.1f\n', new_lp, sqrt(new_var));
    
    % Search again for the new best point
    gp_map(itr+1,:) = gpmax(gp_lp, gp_map(itr,:), [], [], [], [], lb, ub, [], opts);
    
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

function [bestx, bestval] = gpmax(gp, start, varargin)
% Helper to search around a point to find the maximum of a GP
scale = gp.KernelInformation.KernelParameters(1:end-1)';
for i=50:-1:1
    [xval(i,:), fval(i,:)] = fmincon(@(x) -gp.predict(x), start+randn(size(scale)).*scale, varargin{:});
end
[xval(end+1,:), fval(end+1,:)] = fmincon(@(x) -gp.predict(x), start, varargin{:});
[bestval, idxbest] = min(fval);
bestx = xval(idxbest,:);
bestval = -bestval;
end

function newmdl = updateGPRMdl(mdl,Xadd,Yadd,refitHyperParameters)
%UPDATEGPRMDL Helper to 'add' a data point to a GP model without necessarily triggering a slow
%hyperparameter update. Thanks to
%https://www.mathworks.com/matlabcentral/answers/424553-how-to-add-points-to-a-trained-gp
if nargin < 4, refitHyperParameters = true; end
kernelfun = mdl.KernelFunction;
kernelparams = mdl.KernelInformation.KernelParameters;
sigma = mdl.Sigma;
constsig = mdl.ModelParameters.ConstantSigma;
beta = mdl.Beta;
Xall = [mdl.X; Xadd];
Yall = [mdl.Y; Yadd];

if refitHyperParameters
    fitarg = 'exact';
else
    fitarg = 'none';
end

newmdl = fitrgp(Xall, Yall, 'Sigma', sigma, 'Beta', beta, 'KernelParameters', ...
    kernelparams, 'KernelFunction', kernelfun, 'FitMethod', fitarg, 'ConstantSigma', constsig);
end