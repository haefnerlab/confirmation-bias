function [param_set, stim_set, choice_set, trial_set] = GetSubjectDataForFitting(subjectId, kernel_kappa, base_params, datadir)

if ~exist('datadir', 'var'), datadir=fullfile(pwd, '../PublishData'); end
memodir = fullfile(datadir, '../Precomputed');

if ~exist('base_params', 'var') || isempty(base_params)
    % By default, use parameters from qualitative fit that reproduced PK trend
    base_params = Model.newModelParams('model', 'is', 'gamma', 0.1, 'samples', 5, 'updates', 5, ...
        'var_x', 0.1, 'lapse', 0.01, 'temperature', 0.05);
end

% Load subject data. '1' is 'ratio' phase and '2' is 'noise' phase
SubjectDataRatio = LoadAllSubjectData(subjectId, 1, datadir);
SubjectDataNoise = LoadAllSubjectData(subjectId, 2, datadir);

% Get signals from reconstructed stimuli
stim_set{1} = LoadOrRun(@ComputeFrameSignals, {SubjectDataRatio, kernel_kappa}, ...
    fullfile(memodir, ['perFrameSignals-' subjectId '-' num2str(kernel_kappa) '-' SubjectDataRatio.phase '.mat']));
stim_set{2} = LoadOrRun(@ComputeFrameSignals, {SubjectDataNoise, kernel_kappa}, ...
    fullfile(memodir, ['perFrameSignals-' subjectId '-' num2str(kernel_kappa) '-' SubjectDataNoise.phase '.mat']));

% Select sub-threshold trials
[~, trial_set{1}] = GaborThresholdTrials(SubjectDataRatio, 1, .6, .4);
[~, thresh] = GaborAnalysis.getThresholdWindow(SubjectDataNoise, 2, 0.5, 0.7, memodir);
[~, trial_set{2}] = GaborThresholdTrials(SubjectDataNoise, 1, thresh);

% Restrict stim and choices to sub-threshold in each task
stim_set{1} = stim_set{1}(trial_set{1}, :);
stim_set{2} = stim_set{2}(trial_set{2}, :);
choice_set{1} = sign(SubjectDataRatio.choice(trial_set{1})' - 0.5);
choice_set{2} = sign(SubjectDataNoise.choice(trial_set{2})' - 0.5);

% Sub-threshold data in the 'ratio' condition (HSLC) has category information 0.6 since all 6:4 or
% 4:6 trials are included in the analysis. Sub-threshold data in the 'noise' condition (LSHC) has
% category info 0.9 since there was always a single 'pulsed' frame of the opposite category. SI is
% likewise set to 0.6 and 0.9 in each condition, but this is essentially arbitrary since it has no
% effect on the model's behavior for a given log-likelhiood odds; as long as the signal_scale
% mapping from stim to model evidence is fit, the actual value of SI is irrelevant.
param_set(1) = Model.setCategorySensoryInfo(base_params, .6, .9);
param_set(2) = Model.setCategorySensoryInfo(base_params, .9, .6);
end