function [bestfit, train_ll, test_ll, model_names] = IntegratorModelComparisonBySubject(subjectId, phase, sub_threshold, group_by, reference, add_null, datadir)
%% Setup
if nargin < 3, sub_threshold = false; end
if nargin < 4, group_by = 'traintest'; end % or 'model'
if nargin < 5, reference = 'ideal'; end
if nargin < 6, add_null = false; end
if nargin < 7, datadir = fullfile(pwd, '..', 'RawData'); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

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

%% Fit models and get cross-validation scores
prefix = [subjectId '-' num2str(kernel_kappa) '-' SubjectData.phase '-' num2str(sub_threshold)];
[bestfit, train_ll, test_ll, ~, model_names] = IntegratorModelComparison(sigs, choices, 0.98, 100, prefix);
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
if startsWith(group_by, 'm')
    h = bar(est_ll_diff); drawnow;
    errorbar(h(1).XData+h(1).XOffset, est_ll_diff(:,1), sem_ll_diff(:,1), 'ok');
    errorbar(h(2).XData+h(2).XOffset, est_ll_diff(:,2), sem_ll_diff(:,2), 'ok');
    legend({'train', 'test'}, 'Location', 'Northwest');
    set(gca, 'XTick', 1:nModels, 'XTickLabel', model_names);
elseif startsWith(group_by, 't')
    h = bar(est_ll_diff'); drawnow;
    for m=1:length(model_names)
        errorbar(h(m).XData+h(m).XOffset, est_ll_diff(m,:), sem_ll_diff(m,:), 'ok');
    end
    legend(model_names, 'Location', 'Northwest');
    set(gca, 'XTick', 1:2, 'XTickLabel', {'train', 'test'});
end
grid on;
ylabel(['\Delta LL from ' reference ' fit']);
title(['Model Comparison: ' subjectId ' [' upper(SubjectData.phase) ']']);
end