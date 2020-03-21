function [bestfit_params, logpost, loglike, logposterr] = FitModelToData(SubjectData, sensor_noise, model_ype, fit_fields, allow_gamma_neg)

if nargin < 3 || isempty(fit_fields)
    % By default, fit everything that can reasonably be fit
    fit_fields = {'prior_C', 'lapse', 'gamma', 'noise', 'temperature', 'bound', 'updates', 'samples'};
end

%% Process inputs... create base params and subselect fittable fields depending on the model type
switch lower(model_ype)
    case 'is'
        base_params = Model.newModelParams('model', 'is', 'gamma', .1, 'var_x', .1, 'frames', ...
            SubjectData.number_of_images, 'temperature', 0.1, 'updates', 5, 'samples', 5);
        fit_fields = intersect(fit_fields, {'prior_C', 'lapse', 'temperature', 'gamma', 'gamma', ...
            'noise', 'temperature', 'updates', 'samples'});
    case 'vb'
        base_params = Model.newModelParams('model', 'vb-czx', 'gamma', .1, 'var_x', .1, 'frames', ...
            SubjectData.number_of_images, 'temperature', 0.1, 'updates', 5, 'step_size', 0.05);
        fit_fields = intersect(fit_fields, {'prior_C', 'lapse', 'temperature', 'gamma', 'gamma', ...
            'noise', 'temperature', 'updates'});
    case 'itb'
        base_params = Model.newModelParams('model', 'itb', 'gamma', .1, 'var_x', .1, 'frames', ...
            SubjectData.number_of_images, 'temperature', 1, 'noise', 0.05);
        fit_fields = intersect(fit_fields, {'prior_C', 'lapse', 'temperature', 'gamma', 'gamma', ...
            'noise', 'temperature', 'bound'});
    case 'ideal'
        base_params = Model.newModelParams('model', 'ideal', 'var_x', .1, 'frames', ...
            SubjectData.number_of_images, 'temperature', 5);
        fit_fields = intersect(fit_fields, {'prior_C', 'lapse', 'temperature'});
end

%% Convert human data to model data
kernel_kappa = 0.16;
uid = [SubjectData.subjectID '-' num2str(kernel_kappa) '-' SubjectData.phase];
signals = LoadOrRun(@ComputeFrameSignals, {SubjectData, kernel_kappa}, ...
    fullfile('../Precomputed', ['perFrameSignals-' uid '.mat']));

[params_set, stim_set, choice_set] = SubjectDataToModelParams(SubjectData, signals, kernel_kappa, ...
    sensor_noise, base_params);
nonempty = ~cellfun(@isempty, choice_set);
params_set = params_set(nonempty);
stim_set = stim_set(nonempty);
choice_set = choice_set(nonempty);

%% Get priors and bounds
distribs = Fitting.defaultDistributions(fit_fields, true, allow_gamma_neg);
LB  = cellfun(@(f) distribs.(f).lb, fit_fields);
UB  = cellfun(@(f) distribs.(f).ub, fit_fields);
PLB = cellfun(@(f) distribs.(f).plb, fit_fields);
PUB = cellfun(@(f) distribs.(f).pub, fit_fields);

%% Do MAP fitting
bads_options = bads('defaults');
bads_options.UncertaintyHandling = 'yes';

for iRun=10:-1:1
    % Initialize randomly on each iteration
    x0 = PLB + (PUB-PLB).*rand(size(PLB));
    [BESTFIT(iRun,:), NLL(iRun), EXITFLAG(iRun)] = bads(...
        @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(params_set, fit_fields, x), distribs, stim_set, choice_set), ...
        x0, LB, UB, PLB, PUB, [], bads_options);
end

%% Get best model
NLL(EXITFLAG <= 0) = inf;
[~, bestRun] = min(NLL);
bestfit_params = Fitting.setParamsFields(params_set(1), fit_fields, BESTFIT(bestRun, :));

%% Evaluate stats on best model
if nargout >= 2
    [logpost, logpri, ~, logposterr] = Fitting.choiceModelLogProb(bestfit_params, distribs, stim_set, choice_set);
    loglike = logpost - logpri;
    logposterr = sqrt(logposterr);
end

end