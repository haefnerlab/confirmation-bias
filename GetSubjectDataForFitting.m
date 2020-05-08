function [param_set, stim_set, choice_set, trial_set] = GetSubjectDataForFitting(subjectId, kernel_kappa, base_params, datadir)

if ~exist('datadir', 'var'), datadir=fullfile(pwd, '../PublishData'); end
memodir = fullfile(datadir, '../Precomputed');

% These 'phase' flags double as indices here. So stim_set{1} will contain HSLC (ratio phase) data,
% and stim_set{2} will contain LSHC (noise phase) data.
NOISE_PHASE = 2;
RATIO_PHASE = 1;

if ~exist('base_params', 'var') || isempty(base_params)
    % By default, use parameters from qualitative fit that reproduced PK trend
    base_params = Model.newModelParams('model', 'is', 'gamma', 0.1, 'samples', 5, 'updates', 5, ...
        'var_x', 0.1, 'lapse', 0.01, 'temperature', 0.1);
end

SubjectDataRatio = LoadAllSubjectData(subjectId, RATIO_PHASE, datadir);
SubjectDataNoise = LoadAllSubjectData(subjectId, NOISE_PHASE, datadir);

% Get signals from reconstructed stimuli
stim_set{RATIO_PHASE} = LoadOrRun(@ComputeFrameSignals, {SubjectDataRatio, kernel_kappa}, ...
    fullfile(memodir, ['perFrameSignals-' subjectId '-' num2str(kernel_kappa) '-' SubjectDataRatio.phase '.mat']));
stim_set{NOISE_PHASE} = LoadOrRun(@ComputeFrameSignals, {SubjectDataNoise, kernel_kappa}, ...
    fullfile(memodir, ['perFrameSignals-' subjectId '-' num2str(kernel_kappa) '-' SubjectDataNoise.phase '.mat']));

% Select sub-threshold trials
[~, trial_set{RATIO_PHASE}] = GaborThresholdTrials(SubjectDataRatio, RATIO_PHASE, .6, .4);
[~, thresh] = GaborAnalysis.getThresholdWindow(SubjectDataNoise, NOISE_PHASE, 0.5, 0.75, memodir);
[~, trial_set{NOISE_PHASE}] = GaborThresholdTrials(SubjectDataNoise, NOISE_PHASE, thresh);

% Restrict stim and choices to sub-threshold in each task
stim_set{RATIO_PHASE} = stim_set{RATIO_PHASE}(trial_set{RATIO_PHASE}, :);
stim_set{NOISE_PHASE} = stim_set{NOISE_PHASE}(trial_set{NOISE_PHASE}, :);
choice_set{RATIO_PHASE} = sign(SubjectDataRatio.choice(trial_set{RATIO_PHASE})' - 0.5);
choice_set{NOISE_PHASE} = sign(SubjectDataNoise.choice(trial_set{NOISE_PHASE})' - 0.5);

% Sub-threshold data in the 'ratio' condition (HSLC) has category information 0.6 since all 6:4 or
% 4:6 trials are included in the analysis. Sub-threshold data in the 'noise' condition (LSHC) has
% category info 0.9 since there was always a single 'pulsed' frame of the opposite category. SI is
% likewise set to 0.6 and 0.9 in each condition, but this is essentially arbitrary since it has no
% effect on the model's behavior for a given log-likelhiood odds; as long as the signal_scale
% mapping from stim to model evidence is fit, the actual value of SI is irrelevant.
param_set(RATIO_PHASE) = Model.setCategorySensoryInfo(base_params, .6, .9);
param_set(NOISE_PHASE) = Model.setCategorySensoryInfo(base_params, .9, .6);
end