%% Setup

clear;
RATIO_PHASE = 1; % aka HSLC
NOISE_PHASE = 2; % aka LSHC
THRESHOLD = 0.7;
KERNEL_KAPPA = 0.16;
DATADIR = fullfile(pwd, '..', 'PublishData');
MEMODIR = fullfile(pwd, '..', 'Precomputed');

DARK_RED = [204 0 0] / 255;
DARK_BLUE = [32 74 135] / 255;
% fade 50% towards white
LIGHT_RED = DARK_RED * .5 + .5;
LIGHT_BLUE = DARK_BLUE * .5 + .5;

ratioSubjects = arrayfun(@(i) sprintf('BPGTask-subject%02d', i), setdiff(1:15, [1 5 12]), 'UniformOutput', false);
noiseSubjects = arrayfun(@(i) sprintf('BPGTask-subject%02d', i), setdiff(1:15, [1 5]), 'UniformOutput', false);

% Informed is the opposite of naive. Subjects 7 through 9 were authors.
informedSubjects = arrayfun(@(i) sprintf('BPGTask-subject%02d', i), [7 8 9], 'UniformOutput', false);

naiveRatioSubjects = setdiff(ratioSubjects, informedSubjects);
naiveNoiseSubjects = setdiff(noiseSubjects, informedSubjects);
naiveBothSubjects = intersect(naiveRatioSubjects, naiveNoiseSubjects);

% Re-order so that informed subjects are at the end
ratioSubjects = [naiveRatioSubjects informedSubjects];
noiseSubjects = [naiveNoiseSubjects informedSubjects];
bothSubjects = sort([intersect(naiveRatioSubjects, naiveNoiseSubjects) informedSubjects]);
is_naive = ismember(bothSubjects, naiveBothSubjects);

% Compute population-level psycho metrics
N = length(bothSubjects);

for n=N:-1:1
    hyphenSplit = strsplit(bothSubjects{n}, '-');
    shortnames{n} = hyphenSplit{2};
end

% Function handle to evaluate psychometric curve (see @plotPsych from psignifit toolbox)
psychofun = @(psychofit,x) (1-psychofit.Fit(3)-psychofit.Fit(4))*psychofit.options.sigmoidHandle(x,psychofit.Fit(1),psychofit.Fit(2))+psychofit.Fit(4);

%% Before computing PKs and fitting models, we need to (re)compute signal statistics per subject per condition. Slow, but only run once.

for iSubject=1:length(noiseSubjects)
    SubjectData = LoadAllSubjectData(noiseSubjects{iSubject}, NOISE_PHASE, DATADIR);
    sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, KERNEL_KAPPA}, ...
        fullfile(MEMODIR, ['perFrameSignals-' SubjectData.subjectID '-' num2str(KERNEL_KAPPA) '-' SubjectData.phase '.mat']));
end

for iSubject=1:length(ratioSubjects)
    SubjectData = LoadAllSubjectData(ratioSubjects{iSubject}, RATIO_PHASE, DATADIR);
    sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, KERNEL_KAPPA}, ...
        fullfile(MEMODIR, ['perFrameSignals-' SubjectData.subjectID '-' num2str(KERNEL_KAPPA) '-' SubjectData.phase '.mat']));
end

%% Figure 1

% Nothing to do - conceptual figures only created in vector graphics editor

%% Figure 2,3 - data analysis

sem = @(vals) std(vals)/sqrt(N);
for iSubject=length(bothSubjects):-1:1
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, NOISE_PHASE, DATADIR);
    
    % Get PM curve as % correct w.r.t. *unsigned* noise parameter.
    warning off;
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(SubjectData, NOISE_PHASE, 0.5, THRESHOLD, DATADIR);
    warning on;

    % Get %correct at 0 signal and threshold level
    pc_zero_signal(iSubject) = psychofun(psycho_fit, 0);
    threshold(iSubject) = getThreshold(psycho_fit, THRESHOLD);
    
    % Get PM curve w.r.t. frame ratio
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, RATIO_PHASE, DATADIR);
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(SubjectData, RATIO_PHASE, 0.3, 0.7, DATADIR);
    slope = getSlopePC(psycho_fit, 0.5);
    pc_60_40_adjusted(iSubject) = 0.5 + .1 * slope;
end
fprintf('Noise condition mean threshold = %f +/- %f, stdev = %f\n', mean(threshold), sem(threshold), std(threshold));
fprintf('Noise condition mean %%correct at zero signal = %f +/- %f, stdev = %f\n', mean(100*pc_zero_signal), sem(100*pc_zero_signal), std(100*pc_zero_signal));
fprintf('Ratio condition mean %%correct at 6:4 = %f +/- %f, stdev = %f\n', mean(100*pc_60_40_adjusted), sem(100*pc_60_40_adjusted), std(100*pc_60_40_adjusted));

nboot_pvalue = 10000;
GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [RATIO_PHASE NOISE_PHASE], 'exponential', nboot_pvalue, is_naive, DATADIR);
GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [RATIO_PHASE NOISE_PHASE], 'linear', nboot_pvalue, is_naive, DATADIR);

method = 'reg-lr'; % options include 'reg-lr', 'lr', 'exp', or 'lin'. This is just for visualization.
GaborAnalysis.DeltaPK(bothSubjects, [RATIO_PHASE NOISE_PHASE], false, method, DATADIR);

%% Psychometric curves

figure;

% Noise condition
subplot(1,2,1);
hold on;
noises = linspace(-0.8, 0.8, 101);
avg_pm_curve = zeros(size(noises));
for iSubject=1:N
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, NOISE_PHASE, DATADIR);
    [pm_fit, ~, ~, ~] = LoadOrRun(@GaborPsychometric, {SubjectData, -2}, ...
        fullfile(MEMODIR, ['PM-noise-signed-' SubjectData.subjectID '.mat']));
    subject_pm_curve = arrayfun(@(x) psychofun(pm_fit, x), noises);
    avg_pm_curve = avg_pm_curve + (subject_pm_curve - avg_pm_curve) / iSubject;
    plot(noises, subject_pm_curve, 'Color', LIGHT_RED, 'LineWidth', 1);
end
plot(noises, avg_pm_curve, 'Color', DARK_RED, 'LineWidth', 2);
xlabel('noise level (\kappa)');
ylabel('percent chose left');
xlim([-inf inf]);

% Ratio condition
subplot(1,2,2);
hold on;
ratios = linspace(0, 1, 101);
avg_pm_curve = zeros(size(ratios));
for iSubject=1:N
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, RATIO_PHASE, DATADIR);
    [pm_fit, ~, ~, ~] = LoadOrRun(@GaborPsychometric, {SubjectData, RATIO_PHASE}, ...
        fullfile(MEMODIR, ['PM-true_ratio-signed-' SubjectData.subjectID '.mat']));
    subject_pm_curve = arrayfun(@(x) psychofun(pm_fit, x), ratios);
    avg_pm_curve = avg_pm_curve + (subject_pm_curve - avg_pm_curve) / iSubject;
    plot(ratios, subject_pm_curve, 'Color', LIGHT_BLUE, 'LineWidth', 1);
end
plot(ratios, avg_pm_curve, 'Color', DARK_BLUE, 'LineWidth', 2);
xlabel('frame ratio left:right');
ylabel('percent chose left');
xlim([-inf inf]);

%% Supplemental PK cross-validation figure

nFold = 20;
pkLLNoise = zeros(4, N, nFold);
pkLLRatio = zeros(4, N, nFold);
for iSubject=N:-1:1
    %% Noise experiment PK cross-validation
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, NOISE_PHASE, DATADIR);
    
    % Compute per-frame signals
    sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, KERNEL_KAPPA}, ...
        fullfile(MEMODIR, ['perFrameSignals-' SubjectData.subjectID '-' num2str(KERNEL_KAPPA) '-' SubjectData.phase '.mat']));
    
    % Select sub-threshold trials
    [floor, thresh] = GaborAnalysis.getThresholdWindow(SubjectData, NOISE_PHASE, 0.5, THRESHOLD, MEMODIR);
    [~, trials] = GaborThresholdTrials(SubjectData, NOISE_PHASE, thresh, floor);
    sigs = sigs(trials, :);
    choices = SubjectData.choice(trials) == +1;
    
    % Run cross-validated model comparison
    memo_name = ['PK-xValidCompare-' SubjectData.phase '-' SubjectData.subjectID '-' num2str(thresh) '-' num2str(floor) '.mat'];
    [pkLLNoise(1, iSubject, :), pkLLNoise(2, iSubject, :), pkLLNoise(3, iSubject, :), pkLLNoise(4, iSubject, :), ~] = ...
        LoadOrRun(@CustomRegression.xValidatePKModels, {sigs, choices, nFold}, fullfile(MEMODIR, memo_name));
    
    %% Ratio experiment PK cross-validation
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, RATIO_PHASE, DATADIR);
    
    % Compute per-frame signals
    sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, KERNEL_KAPPA}, ...
        fullfile(MEMODIR, ['perFrameSignals-' SubjectData.subjectID '-' num2str(KERNEL_KAPPA) '-' SubjectData.phase '.mat']));
    
    % Select sub-threshold trials
    [~, trials] = GaborThresholdTrials(SubjectData, RATIO_PHASE, .6, .4);
    sigs = sigs(trials, :);
    choices = SubjectData.choice(trials) == +1;
    
    % Run cross-validated model comparison
    memo_name = ['PK-xValidCompare-' SubjectData.phase '-' SubjectData.subjectID '-0.6-0.4.mat'];
    [pkLLRatio(1, iSubject, :), pkLLRatio(2, iSubject, :), pkLLRatio(3, iSubject, :), pkLLRatio(4, iSubject, :), ~] = ...
        LoadOrRun(@CustomRegression.xValidatePKModels, {sigs, choices, nFold}, fullfile(MEMODIR, memo_name));
end

% Note: CustomRegression.xValidatePKModels gives paired results. That is, each of the 4 models is
% tested on the exact same data fold. The relevant statistic is therefore the *average of
% differences in LL* rather than the *difference in averages*. Size of both 'pkLLNoise' and
% 'pkLLRatio' is [#models x #subjects x #folds]
xValSubjects = shortnames;
models = {'LR', 'reg. LR', 'Exp.', 'Lin.'};

% Sum LL across subjects to get LL for the whole dataset
pkLLNoise(:, end+1, :) = sum(pkLLNoise, 2);
pkLLRatio(:, end+1, :) = sum(pkLLRatio, 2);
xValSubjects{end+1} = 'Combined';

% Compute differences w.r.t. unregularized logistic regression (essentially glmfit) as baseline
baselineModel = 1;
pkLLNoiseDiffs = pkLLNoise - pkLLNoise(baselineModel, :, :);
pkLLRatioDiffs = pkLLRatio - pkLLRatio(baselineModel, :, :);

% Compute mean and standard error across folds
meanLLNoise = mean(pkLLNoiseDiffs, 3);
semLLNoise = std(pkLLNoiseDiffs, [], 3) ./ sqrt(nFold);
meanLLRatio = mean(pkLLRatioDiffs, 3);
semLLRatio = std(pkLLRatioDiffs, [], 3) ./ sqrt(nFold);

figure;
subplot(2,1,1);
hold on;
p1 = bar(meanLLNoise');
drawnow;
for iSubject=1:length(xValSubjects)
    for iModel=1:length(models)
        x = p1(iModel).XData(iSubject) + p1(iModel).XOffset;
        errorbar(x, meanLLNoise(iModel, iSubject), semLLNoise(iModel,iSubject), 'Color', 'k');
    end
end
legend({'Unregularized LR', 'Regularized LR', 'Exponential', 'Linear'}, 'location', 'best');
ylabel('Relative log-likelihood');
title('Noise Condition');
set(gca, 'XTick', 1:length(xValSubjects), 'XTickLabel', xValSubjects);
xtickangle(45);

subplot(2,1,2);
hold on;
p1 = bar(meanLLRatio');
drawnow;
for iSubject=1:length(xValSubjects)
    for iModel=1:length(models)
        x = p1(iModel).XData(iSubject) + p1(iModel).XOffset;
        errorbar(x, meanLLRatio(iModel, iSubject), semLLRatio(iModel,iSubject), 'Color', 'k');
    end
end
legend({'Unregularized LR', 'Regularized LR', 'Exponential', 'Linear'}, 'location', 'best');
ylabel('Relative log-likelihood');
title('Ratio Condition');
set(gca, 'XTick', 1:length(xValSubjects), 'XTickLabel', xValSubjects);
xtickangle(45);

%% Figure 4 - simulation of hierarchical ITB model

fig4 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

beta_range = [-0.5 .15];
pk_hprs = [1 0 0];

% Case 1: no bound, no leak (aka ideal observer)
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'temperature', 0.05, 'bound', inf, 'gamma', 0, 'noise', 0.35);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 1, parula);
figureToPanel(pk_fig, fig4, 5, 3, 2);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 3, cmap);

% Case 2: bound but no leak (Kiani et al 2008)
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'temperature', 0.05, 'bound', 1.2, 'gamma', 0, 'noise', 0.35);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 4, parula);
figureToPanel(pk_fig, fig4, 5, 3, 5);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 6, cmap);

% Case 3: leak but no bound
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'temperature', 0.05, 'bound', inf, 'gamma', .1, 'noise', 0.35);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 7, parula);
figureToPanel(pk_fig, fig4, 5, 3, 8);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 9, cmap);

% Case 4: both leak and bound
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'temperature', 0.05, 'bound', 1.2, 'gamma', .1, 'noise', 0.35);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 10, parula);
figureToPanel(pk_fig, fig4, 5, 3, 11);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 12, cmap);

% Case 5: like previous but now gamma changes with category info
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'temperature', 0.05, 'bound', 1.2, 'gammafun', @(ci,si) (1-ci), 'noise', 0.35);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 13, parula);
figureToPanel(pk_fig, fig4, 5, 3, 14);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 15, cmap);

%% Figure 5 - IS and VB simulation results

fig5 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% Ideal observer performance
params = Model.newModelParams('model', 'ideal', 'trials', 10000);
disp('Loading/Running ideal observer');
[~, ideal_fig] = Model.plotCategorySensorySpace(ps, ps, params);
figureToPanel(ideal_fig, fig5, 2, 4, 5, parula);

% Sampling model with gamma = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)
disp('Loading/Running sampling model, getting threshold points for PKs');
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig5, 2, 4, 2, parula);
figureToPanel(pk_fig, fig5, 2, 4, 3);
disp('Loading/Running sampling model, gettings slopes over CS-Space');
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, fig5, 2, 4, 4, cmap);

% Variational model with gamma = 0.1
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)
disp('Loading/Running variational model, getting threshold points for PKs');
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig5, 2, 4, 6, parula);
figureToPanel(pk_fig, fig5, 2, 4, 7);
disp('Loading/Running variational model, gettings slopes over CS-Space');
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, fig5, 2, 4, 8, cmap);

%% Figure 6 - model optimality figure varying 'gamma' parameter

fig6 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% >> Uncomment for sampling model <<
lo_gamma = 0.1;
hi_gamma = 0.5;
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', hi_gamma, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'samples', 5);

% >> Uncomment for variational model <<
% lo_gamma = 0.1;
% hi_gamma = 0.5;
% params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', hi_gamma, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'step_size', 0.05);

% First panel: percent correct with high gamma (low gamma was done in Fig 3)
[correct_hi, cs_fig] = Model.plotCategorySensorySpace(ps, ps, params);
figureToPanel(cs_fig, fig6, 2, 3, 1, parula);
title(sprintf('Performance (\\gamma=%.1f)', hi_gamma));

% Second panel: difference in % correct
params_lo = setfield(params, 'gamma', lo_gamma);
[~, diff_fig] = Model.plotDeltaPerformance(ps, ps, params, params_lo);
figureToPanel(diff_fig, fig6, 2, 3, 2, parula);
title(sprintf('\\DeltaPerformance (\\gamma=%.1f vs. \\gamma=%.1f)', hi_gamma, lo_gamma));

% Third panel: optimal gamma across C-S space
[~, opt_correct_fig, opt_gamma_fig] = Model.plotCategorySensorySpace(ps, ps, params, {'gamma'}, 21);
figureToPanel(opt_gamma_fig, fig6, 2, 3, 3, cool);
title('Optimized value of \gamma');

% Fourth panel: percent correct with optimal gamma
figureToPanel(opt_correct_fig, fig6, 2, 3, 4, parula);
title('Performance with optimal \gamma*');

% Fifth panel: diff in % correct between optimized model and ideal observer
params_ideal = setfield(params, 'model', 'ideal');
[~, diff_fig] = Model.plotDeltaPerformance(ps, ps, params_ideal, params, {}, 0, {'gamma'}, 21);
figureToPanel(diff_fig, fig6, 2, 3, 5, parula);
title('\DeltaPerformance (ideal vs. \gamma*)');

% Sixth panel: slopes with optimal gamma
beta_range = [-.1 .1];
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, [], {'gamma'}, 21);
figureToPanel(fig, fig6, 2, 3, 6, cmap);
title('Temporal Weight Slope \beta with \gamma*');

%% Supplemental figure model results (showing effect of gamma parameter)

figModelSupp = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% --- SAMPLING ---

% Replicate sampling model results above (gamma = 0)
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-0.4 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 1, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 2);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 3, cmap);

% Same with gamma = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.43 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 4, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 5);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 6, cmap);

% Same with gamma = 0.2
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.2, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.4 .23]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 7, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 8);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 9, cmap);

% --- VB-CZX ---

% Replicate vb model results above (gamma = 0)
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-0.5 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 10, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 11);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 12, cmap);

% Same with gamma = 0.1
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 13, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 14);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 15, cmap);

% Same with gamma = 0.2
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.2, 'temperature', 0.05, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.26 .2]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 16, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 17);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 18, cmap);

%% Pre-run MCMC sampling of model parameters

% NOTE: this will take a very very long time. In practice, we run this on a cluster, parallelized
% across chains, for weeks to actually generate MCMC chains, and subsequent sections are designed to
% load and plot "whatever is done so far" (setting nSamplesPerChain to 0 will load existing samples
% but not generate new ones). This section is included here primarily as documentation.

% ITB and IS are ground-truth models (see @GetGroundTruthSimData)
subjectsToFit = [{'ITB', 'IS'} bothSubjects];
phasesToFit = {'lshc', 'hslc', 'both'};
nSamplesPerChain = 5e5;
nChains = 12;
for iSub=1:length(subjectsToFit)
    for iPhz=1:length(phases)
        GetITBPosteriorSamples(subjectsToFit{iSub}, phases{iPhz}, nSamplesPerChain, nChains, DATADIR, MEMODIR);
    end
end

%% Supplemental figure: all inferred parameter values

subjectsToFit = [{'ITB', 'IS'} bothSubjects];
fields = {'prior_C', 'lapse', 'temperature', 'signal_scale', 'neggamma', 'bound', 'noise'};
colors = [lines(4); 0 .5 0; .4 0 .4; .4 0 .4];
[fig_scatter, fig_histogram, fig_para_box] = ITBParameterPlot(subjectsToFit, {'lshc', 'hslc'}, ...
    fields, [0 1 1 1 0 1 1], colors, 'oooo^s*', DATADIR, MEMODIR);

% Keep the 'box' style plot; close the other two
close(fig_scatter);
close(fig_histogram);

%% Supplemental bar plots of 'Beta Explained' for all subjects and models

[tmp, fig_beta_bar_gt] = BetaExplained({'ITB', 'IS'}, 'gbn', DATADIR, MEMODIR);
close(tmp);
[tmp, fig_beta_bar_subjects, betaDiagnostics, phase_names] = BetaExplained(bothSubjects, 'gbn', DATADIR, MEMODIR);
close(tmp);

% Print out convergence info for ablation tests
for iPhz=1:length(phase_names)
    maxRHat = max(cellfun(@(info) max(info.RHat), betaDiagnostics(:,iPhz)));
    minRHat = min(cellfun(@(info) min(info.RHat), betaDiagnostics(:,iPhz)));
    fprintf('Phase %s beta rhat in [%.3f %.3f]\n', upper(phase_names{iPhz}), minRHat, maxRHat);
end

%% Model fitting analysis on subjects

fig_fit_results = figure;

% Panels 1 and 2 along top third: scatter plots of the gamma and bound parameters
fields = {'neggamma', 'bound'};
colors = [0 .5 0; .4 0 .4];
[fig_param_scatter, fig_histogram, fig_box] = ITBParameterPlot(bothSubjects, {'lshc', 'hslc'}, ...
    fields, [0 1 1 0 1 1 1], colors, '^s', DATADIR, MEMODIR);
close(fig_histogram);
close(fig_box);
figureToPanel(fig_param_scatter, fig_fit_results, 2, 2, 1);

fig_beta_scatter = BetaExplained(bothSubjects, '', DATADIR, MEMODIR);
for iax=1:length(fig_beta_scatter.Children)
    if isaxes(fig_beta_scatter.Children(iax))
        % Delete all but 'full model' series. Children(1:2) are subplots, each of which has
        % Children(1:3) that are 3 errorbar series
        delete(fig_beta_scatter.Children(iax).Children(1:2));
    end
end
% Panels 3 and 4 along top third: 'Beta Explained' of full model for each condition.
figureToPanel(fig_beta_scatter, fig_fit_results, 2, 2, 2);

% Bottom 2 sets of panels: another copy of fig_beta_scatter but isolating the gamma and bound terms
fig_beta_scatter = BetaExplained(bothSubjects, '', DATADIR, MEMODIR);
for iax=1:length(fig_beta_scatter.Children)
    if isaxes(fig_beta_scatter.Children(iax))
        % Delete full-model series
        delete(fig_beta_scatter.Children(iax).Children(3));
    end
end
% Panels 3 and 4 along top third: 'Beta Explained' of full model for each condition.
figureToPanel(fig_beta_scatter, fig_fit_results, 2, 1, 2);

%% Helper function for figure layout

function ax_copy = figureToPanel(figSource, figDest, subM, subN, subI, customcmap)
margin = 0.02;
widths = (1 - (subN + 1) * margin) / subN;
heights = (1 - (subM + 1) * margin) / subM;
[col, row] = ind2sub([subN, subM], subI);
left = margin + (col-1) * (widths + margin);
top = margin + (row-1) * (heights + margin);
bottom = 1 - top - heights;
figure(figSource);
for iax=1:length(figSource.Children)
    ax = figSource.Children(iax);
    if ~isaxes(ax), continue; end
    ax_copy = copyobj(ax, figDest);
    new_l = left + widths*ax.Position(1);
    new_b = bottom + heights*ax.Position(2);
    new_w = widths*ax.Position(3);
    new_h = heights*ax.Position(4);
    ax_copy.Position = [new_l new_b new_w new_h];
    ax_copy.Parent = figDest;
    if exist('customcmap', 'var')
        colormap(ax_copy, customcmap);
        colorbar('peer', ax_copy);
    end
end
close(figSource);
drawnow;
end

function tf = isaxes(ax)
%https://www.mathworks.com/matlabcentral/answers/300880-what-is-best-practice-to-determine-if-input-is-a-figure-or-axes-handle
try
    tf = strcmp(get(ax, 'type'), 'axes');
catch
    tf = false;
end
end