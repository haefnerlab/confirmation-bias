function [samples_by_model, model_names, model_fields, params_samples] = IntegratorSampleBySubject(subjectId, phase, n_samples, sub_threshold, datadir, cache)
%% Setup
if nargin < 3, n_samples = 10000; end
if nargin < 4, sub_threshold = false; end
if nargin < 5, datadir = fullfile(pwd, '..', 'RawData'); end
if nargin < 6, cache = true; end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

%% Define models
model_fields = {
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'bound'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'gamma'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'bound', 'gamma'}
};

model_names = {'ideal', 'itb', 'gamma', 'itb_gamma'};

%% Load subject data
SubjectData = LoadAllSubjectData(subjectId, phase, datadir);
kernel_kappa = 0.16;

sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, kernel_kappa}, ...
    fullfile(memodir, ['perFrameSignals-' SubjectData.subjectID '-' num2str(kernel_kappa) '-' SubjectData.phase '.mat']));
choices = SubjectData.choice(:) == +1;

if sub_threshold
    if phase == 1
        thresh = 0.6;
        floor = 0.4;
    elseif phase == 2
        [floor, thresh] = GaborAnalysis.getThresholdWindow(SubjectData, phase, 0.5, 0.7, memodir);
    else
        error('not implemented');
    end
    [~, trials] = GaborThresholdTrials(SubjectData, phase, thresh, floor);
    sigs = sigs(trials, :);
    choices = choices(trials, :);
end

%% Run sampling
phaseName = SubjectData.phase;
parfor iModel=1:length(model_names)
    if cache
        uid = ['integrator-' model_names{iModel} '-samples-' subjectId '-' num2str(kernel_kappa) '-' phaseName '-' num2str(sub_threshold) '-ns=' num2str(n_samples)];
        [samples_by_model{iModel}, ~, base_params] = LoadOrRun(@Fitting.sampleIntegratorModel, ...
            {sigs, choices, n_samples, model_fields{iModel}, true, [], 'burnin', 500, 'thin', 20}, ...
            fullfile(memodir, [uid '.mat']));
    else
        [samples_by_model{iModel}, ~, base_params] = Fitting.sampleIntegratorModel(sigs, choices, ...
            n_samples, model_fields{iModel}, true, [], 'burnin', 500, 'thin', 20);
    end
    params_samples{iModel} = repmat(base_params, n_samples, 1);
    for s=1:n_samples
        params_samples{iModel}(s) = Fitting.setParamsFields(base_params, model_fields{iModel}, samples_by_model{iModel}(s,:));
    end
end

end