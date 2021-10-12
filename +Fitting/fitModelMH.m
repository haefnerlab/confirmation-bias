function [optim_results, samples, sample_scores, metadata] = fitModelMH(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELMH Strategy for fitting inference models by MH sampling followed by selecting the
%best sample as an approximation to the MAP/MLE. Stochastic models are handled as a special case.

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
chkpt = fullfile('fit-checkpoints', sprintf('%x', input_id));
ibs_repeats = 5;

nPara = length(base_params);

%% Run MH sampling and handle stochastic vs deterministic cases
if Model.isStochastic(base_params)
    % Note that stochastic models with finite repetitions will result in posteriors that are wider
    % than the true posterior because each LH evaluation is corrupted by noise
    disp('fitModelMH: draw samples (stochastic case)');
    [samples, ~, sample_scores] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, ibs_repeats, 50, 10, chkpt);

    gp_sigma = sqrt(mean(sample_scores.variance));
    eval_fun = @Fitting.choiceModelLogProbIBS;
    eval_args = {[], ibs_repeats};
    final_eval_args = {[], round(10*sqrt(nTrials))};
else
    disp('fitModelMH: draw samples (deterministic case)');
    [samples, ~, sample_scores] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 0, 50, 10, chkpt);
    
    gp_sigma = 0.01;
    eval_fun = @Fitting.choiceModelLogProb;
    eval_args = {1};
    final_eval_args = {1};
end

%% Find best sample based on existing (rough) scores
disp('fitModelMH: search samples');

[uSamples, ~, idxExpand] = unique(samples, 'rows');
for iSamp=size(uSamples, 1):-1:1
    combo_ll(iSamp) = mean(sample_scores.loglike(idxExpand == iSamp));
    combo_lp(iSamp) = mean(sample_scores.logpdf(idxExpand == iSamp));
end

[~, best_ll_idx] = max(combo_ll);
optim_results{1} = Fitting.setParamsFields(base_params, fields, uSamples(best_ll_idx, :));
for iP=1:nPara
    optim_results{1}(iP).fit_fields = fields;
    optim_results{1}(iP).fit_method = 'mh-mle';
end

[~, best_lp_idx] = max(combo_lp);
optim_results{2} = Fitting.setParamsFields(base_params, fields, uSamples(best_lp_idx, :));
for iP=1:nPara
    optim_results{2}(iP).fit_fields = fields;
    optim_results{2}(iP).fit_method = 'mh-map';
end

%% Search LL landscape with GP
disp('fitModelMH: GP MLE');

% Set initial length scales as 1/20th of the range of plausible bounds for each parameter
init_scale = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/20, fields)';
% Set convergence tolerance as 1/1000th of the plausible range
tolerance = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/1000, fields)';

% Fit Gaussian process to all log likelihood evaluations
gp_ll = fitrgp(samples, sample_scores.loglike, 'KernelFunction', 'ardsquaredexponential', ...
    'KernelParameters', [init_scale std(combo_ll)], 'Sigma', gp_sigma);

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
lb = cellfun(@(f) distribs.(f).lb, fields);
ub = cellfun(@(f) distribs.(f).ub, fields);
gp_mle(1,:) = gpmax(gp_ll, uSamples(best_ll_idx, :), [], [], [], [], lb, ub, [], opts);

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_ll.predict(gp_mle(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_mle=%.1f+/-%.1f', itr, max(abs(delta)), ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [~, new_ll, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_mle(itr,:)), distribs, signals, choices, eval_args{:});
    refit = mod(itr,20)==0;
    gp_ll = updateGPRMdl(gp_ll, gp_mle(itr,:), new_ll, refit);
    
    fprintf('\tactual_mle=%.1f+/-%.1f\n', new_ll, sqrt(new_var));
    
    % Search again for the new best point
    gp_mle(itr+1,:) = gpmax(gp_ll, gp_mle(itr,:), [], [], [], [], lb, ub, [], opts);
    
    delta = gp_mle(itr+1,:) - gp_mle(itr,:);
    itr = itr+1;
end

optim_results{3} = Fitting.setParamsFields(base_params, fields, gp_mle(end,:));
for iP=1:nPara
    optim_results{3}(iP).fit_fields = fields;
    optim_results{3}(iP).fit_method = 'mh-gp-mle-search';
end

%% Repeat the above for the MAP
disp('fitModelMH: GP MAP');

% Fit Gaussian process to all log likelihood evaluations
gp_lp = fitrgp(samples, sample_scores.logpdf, 'KernelFunction', 'ardsquaredexponential', ...
    'KernelParameters', [init_scale std(combo_lp)], 'Sigma', gp_sigma);

% Search the GP fit for a better maximum
opts = optimoptions('fmincon', 'display', 'none');
gp_map(1,:) = gpmax(gp_lp, uSamples(best_lp_idx, :), [], [], [], [], lb, ub, [], opts);

itr = 1;
delta = inf;
while ~all(abs(delta) < tolerance)
    [ypred, ypred_sd] = gp_lp.predict(gp_map(itr,:));
    fprintf('Iteration %d\tdelta=%.1e\test_logp=%.1f+/-%.1f', itr, max(abs(delta)), ypred, ypred_sd);
    
    % Evaluate the current max point and update the GP
    [new_lp, ~, new_var] = eval_fun(Fitting.setParamsFields(base_params, fields, gp_map(itr,:)), distribs, signals, choices, eval_args{:});
    refit = mod(itr,20)==0;
    gp_lp = updateGPRMdl(gp_lp, gp_map(itr,:), new_lp, refit);
    
    fprintf('\tactual_logp=%.1f+/-%.1f\n', new_lp, sqrt(new_var));
    
    % Search again for the new best point
    gp_map(itr+1,:) = gpmax(gp_lp, gp_map(itr,:), [], [], [], [], lb, ub, [], opts);
    
    delta = gp_map(itr+1,:) - gp_map(itr,:);
    itr = itr+1;
end

optim_results{4} = Fitting.setParamsFields(base_params, fields, gp_map(end,:));
for iP=1:nPara
    optim_results{4}(iP).fit_fields = fields;
    optim_results{4}(iP).fit_method = 'mh-gp-map-search';
end

%% Re-run and store final evaluations for each model
for iFit=1:length(optim_results)
    fprintf('fitModelMH: final evals [%s]\n', optim_results{iFit}(1).fit_method);
    [~, ll, ll_var] = eval_fun(optim_results{iFit}, distribs, signals, choices, final_eval_args{:});
    for iP=1:nPara
        optim_results{iFit}(iP).ll = ll;
        optim_results{iFit}(iP).ll_var = ll_var;
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