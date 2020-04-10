function [aic, aic_err, model_info, sampleses] = ModelComparison(base_params, signals, choices, fit_scale, prefix, model_names)

% Note: use 'base_params' to set generative model parameters like CI, SI, etc

%% Specify models
% Each struct in the array is a model to fit.
model_info = struct(...
    'name', {'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split', 'ideal'}, ...
    'type', {'is', 'vb-czx', 'itb', 'itb', 'itb', 'itb', 'ideal'}, ...
    'fields', {{'prior_C', 'log_lapse', 'gamma', 'samples'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'step_size', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma_1', 'gamma_2', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma_1', 'gamma_2', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse'}}, ...
    'gammaneg', {0 0 0 1 0 1 0});

if length(base_params) == 1
    % When fitting only a single condition, remove the 'split' models
    model_info = model_info(cellfun(@(nm) ~contains(nm, 'split'), {model_info.name}));
end

% If supplied, restrict to the ones requested
if nargin >= 5
    names = {model_info.name};
    keep = cellfun(@(nm) any(strcmpi(nm, model_names)), names);
    model_info = model_info(keep);
end

%% Fit each model
use_cache =  nargin >= 4 && ~isempty(prefix);
parfor iModel=1:length(model_info)
    this_params = base_params;
    fields = model_info(iModel).fields;
    if fit_scale
        fields = [fields {'log_signal_scale'}];
    end
    for iP=1:length(this_params)
        % params may be a struct array
        this_params(iP).model = model_info(iModel).type;
    end
    distribs = Fitting.defaultDistributions(fields, false, model_info(iModel).gammaneg);
    % Randomly init all params. This is primarily so that Model.isStochastic() returns true when
    % 'noise' is a parameter even if base_params.noise=0. But it also just seems like good practice.
    this_params = Fitting.setParamsFields(this_params, fields, cellfun(@(f) distribs.(f).priorrnd(1), fields));
    if use_cache
        [optim_results, sampleses{iModel}, ~, ~] = LoadOrRun(@Fitting.fitModelMH, ...
            {this_params, signals, choices, distribs, struct('prefix', prefix)}, ...
            fullfile('../Precomputed', ['mhfit-' prefix '-' model_info(iModel).name '.mat']));
    else
        [optim_results, sampleses{iModel}, ~, ~] = Fitting.fitModelMH(this_params, signals, choices, distribs, struct('prefix', prefix));
    end
    
    % TODO - cache or do this inside of fitModelMH
    if iscell(choices)
        nTrials = sum(cellfun(@length, choices));
    else
        nTrials = length(choices);
    end
    [~, mle(iModel), ll_var(iModel)] = Fitting.choiceModelLogProbIBS(optim_results.gp_mle_params, distribs, signals, choices, [], round(sqrt(nTrials)));
    npara(iModel) = length(fields);
end

%% Compute AIC
aic = aicbic(mle, npara);
aic_err = sqrt(ll_var);

end