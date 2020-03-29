function [bestfit, train_ll, test_ll, model_names] = IntegratorModelComparisonBySubject(subjectId, phase, sub_threshold, reference, add_null, datadir)
%% Setup
if nargin < 3, sub_threshold = false; end
if nargin < 4, reference = 'ideal'; end
if nargin < 5, add_null = false; end
if nargin < 6, datadir = fullfile(pwd, '..', 'RawData'); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

%% Load subject data
SubjectData = LoadAllSubjectData(subjectId, phase, '../PublishData');
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

%% Fit models and get cross-validation scores
prefix = [subjectId '-' num2str(kernel_kappa) '-' SubjectData.phase];
[bestfit, train_ll, test_ll, model_names] = IntegratorModelComparison(sigs, choices, 50, prefix, 989895774);
nModels = length(model_names);

if add_null
    % For reference, compute the log likelihood under the null model
    ll_null = -log(2);
    nModels = nModels+1;
    model_names = [model_names {'null'}];
    train_ll(nModels, :) = ll_null;
    test_ll(nModels, :) = ll_null;
end

% Select which model to use as baseline
ref_model = strcmpi(model_names, reference);

%% Compute LL _differences_ across reruns

% Difference on training sets
ll_diff = train_ll - train_ll(ref_model,:);
est_ll_diff(:,1) = mean(ll_diff, 2);
sem_ll_diff(:,1) = std(ll_diff, [], 2) ./ sqrt(size(ll_diff, 2));

% Difference on validation sets
ll_diff = test_ll - test_ll(ref_model,:);
est_ll_diff(:,2) = mean(ll_diff, 2);
sem_ll_diff(:,2) = std(ll_diff, [], 2) ./ sqrt(size(ll_diff, 2));

%% Plot result
hold on;
bar(1:nModels, est_ll_diff);
errorbar((1:nModels)-.15, est_ll_diff(:,1), sem_ll_diff(:,1), 'ok');
errorbar((1:nModels)+.15, est_ll_diff(:,2), sem_ll_diff(:,2), 'ok');
set(gca, 'XTick', 1:nModels, 'XTickLabel', model_names);
grid on;
ylabel(['\Delta LL relative to ' reference]);
title(['Model Comparison: ' subjectId ' [' upper(SubjectData.phase) ']']);
end