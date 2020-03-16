%% Setup

clear;
RATIO_PHASE = 1; % aka HSLC
NOISE_PHASE = 2; % aka LSHC
THRESHOLD = 0.7;
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
bothSubjects = [intersect(naiveRatioSubjects, naiveNoiseSubjects) informedSubjects];
is_naive = ismember(bothSubjects, naiveBothSubjects);

% Compute population-level psycho metrics
N = length(bothSubjects);

for n=N:-1:1
    hyphenSplit = strsplit(bothSubjects{n}, '-');
    shortnames{n} = hyphenSplit{2};
end

%% Figure 1

% Nothing to do - conceptual figures only created in vector graphics editor

%% Figure 2,3 - data analysis

sem = @(vals) std(vals)/sqrt(N);
for iSubject=length(bothSubjects):-1:1
    warning off;
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, NOISE_PHASE, 0.5, THRESHOLD, DATADIR);
    warning on;
    % See @plotPsych, e.g.
    pc_zero_signal(iSubject) = (1-psycho_fit.Fit(3)-psycho_fit.Fit(4))*psycho_fit.options.sigmoidHandle(0,psycho_fit.Fit(1),psycho_fit.Fit(2))+psycho_fit.Fit(4);
    threshold(iSubject) = getThreshold(psycho_fit, 0.75);
    
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, RATIO_PHASE, 0.3, 0.7, DATADIR);
    slope = getSlopePC(psycho_fit, 0.5);
    pc_60_40_adjusted(iSubject) = 0.5 + .1 * slope;
end
fprintf('Noise condition mean threshold = %f +/- %f, stdev = %f\n', mean(threshold), sem(threshold), std(threshold));
fprintf('Noise condition mean %%correct at zero signal = %f +/- %f, stdev = %f\n', mean(100*pc_zero_signal), sem(100*pc_zero_signal), std(100*pc_zero_signal));
fprintf('Ratio condition mean %%correct at 6:4 = %f +/- %f, stdev = %f\n', mean(100*pc_60_40_adjusted), sem(100*pc_60_40_adjusted), std(100*pc_60_40_adjusted));

nboot_pvalue = 10000;
GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [1 2], 'exponential', nboot_pvalue, is_naive, DATADIR);
GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [1 2], 'linear', nboot_pvalue, is_naive, DATADIR);

method = 'reg-lr'; % options include 'reg-lr', 'lr', 'exp', or 'lin'
GaborAnalysis.DeltaPK(bothSubjects, [1 2], false, method, DATADIR);

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
        fullfile(MEMODIR, ['PM-noise-signed-' bothSubjects{iSubject} '.mat']));
    % Next line copied from @plotPsych in the psignifit toolbox
    subject_pm_curve = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), noises)+pm_fit.Fit(4);
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
        fullfile(MEMODIR, ['PM-true_ratio-signed-' bothSubjects{iSubject} '.mat']));
    % Next line copied from @plotPsych in the psignifit toolbox
    subject_pm_curve = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
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
    % Noise experiment PK cross-validation
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, NOISE_PHASE, DATADIR);
    [floor, thresh] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, NOISE_PHASE, 0.5, THRESHOLD);
    SubjectDataThresh = GaborThresholdTrials(SubjectData, NOISE_PHASE, thresh, floor);
    sigs = SubjectDataThresh.ideal_frame_signals;
    resps = SubjectDataThresh.choice == +1;
    memo_name = ['PK-xValidCompare-noise-' bothSubjects{iSubject} '-' num2str(thresh) '-' num2str(floor) '.mat'];
    [pkLLNoise(1, iSubject, :), pkLLNoise(2, iSubject, :), pkLLNoise(3, iSubject, :), pkLLNoise(4, iSubject, :), ~] = ...
        LoadOrRun(@CustomRegression.xValidatePKModels, {sigs, resps, nFold}, fullfile(MEMODIR, memo_name));
    
    % Ratio experiment PK cross-validation
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, RATIO_PHASE, DATADIR);
    SubjectDataThresh = GaborThresholdTrials(SubjectData, RATIO_PHASE, .6, .4);
    sigs = SubjectDataThresh.ideal_frame_signals;
    resps = SubjectDataThresh.choice == +1;
    disp(bothSubjects{iSubject});
    disp(size(sigs));
    disp(size(resps));
    memo_name = ['PK-xValidCompare-true_ratio-' bothSubjects{iSubject} '-0.6-0.4.mat'];
    [pkLLRatio(1, iSubject, :), pkLLRatio(2, iSubject, :), pkLLRatio(3, iSubject, :), pkLLRatio(4, iSubject, :), ~] = ...
        LoadOrRun(@CustomRegression.xValidatePKModels, {sigs, resps, nFold}, fullfile(MEMODIR, memo_name));
end

% Note: CustomRegression.xValidatePKModels gives paired results. That is, each of the 4 models is
% tested on the exact same data fold. The relevant statistic is therefore the *average of
% differences in LL* rather than the *difference in averages*.

% Compute differences w.r.t. linear model as the baseline
baselineModel = 4;
pkLLNoiseDiffs = pkLLNoise - pkLLNoise(baselineModel, :, :);
pkLLRatioDiffs = pkLLRatio - pkLLRatio(baselineModel, :, :);

% Sum across subjects for combined result
pkLLNoiseDiffsAll = squeeze(sum(pkLLNoiseDiffs, 2));
pkLLRatioDiffsAll = squeeze(sum(pkLLNoiseDiffs, 2));

% Sort across 'folds' so it's easy to get confidence intervals
pkLLNoiseDiffs = sort(pkLLNoiseDiffs, 3);
pkLLRatioDiffs = sort(pkLLRatioDiffs, 3);
pkLLNoiseDiffsAll = sort(pkLLNoiseDiffsAll, 2);
pkLLRatioDiffsAll = sort(pkLLRatioDiffsAll, 2);

% Get mean and 50% confidence intervals across folds for each statistic
meanDiffLLNoise = mean(pkLLNoiseDiffs, 3);
loLLNoise = pkLLNoiseDiffs(:, :, round(.25 * nFold));
hiLLNoise = pkLLNoiseDiffs(:, :, round(.75 * nFold));

meanDiffLLRatio = mean(pkLLRatioDiffs, 3);
loLLRatio = pkLLRatioDiffs(:, :, round(.25 * nFold));
hiLLRatio = pkLLRatioDiffs(:, :, round(.75 * nFold));

meanDiffLLNoiseAll = mean(pkLLNoiseDiffsAll, 2);
loLLNoiseAll = pkLLNoiseDiffsAll(:, round(.25 * nFold));
hiLLNoiseAll = pkLLNoiseDiffsAll(:, round(.75 * nFold));

meanDiffLLRatioAll = mean(pkLLRatioDiffsAll, 2);
loLLRatioAll = pkLLRatioDiffsAll(:, round(.25 * nFold));
hiLLRatioAll = pkLLRatioDiffsAll(:, round(.75 * nFold));

figure;
subplot(2,1,1);
hold on;
p1 = bar([meanDiffLLNoise'; meanDiffLLNoiseAll']);
drawnow;
for iSubject=1:N
    for iType=1:4
        c = p1(iType).FaceColor;
        x = p1(iType).XData(iSubject) + p1(iType).XOffset;
        m = meanDiffLLNoise(iType, iSubject);
        errorbar(x, m, m-loLLNoise(iType, iSubject), hiLLNoise(iType, iSubject)-m, 'Color', 'k');
    end
end
for iType=1:4
    c = p1(iType).FaceColor;
    x = N + 1 + p1(iType).XOffset;
    m = meanDiffLLNoiseAll(iType);
    errorbar(x, m, m-loLLNoiseAll(iType), hiLLNoiseAll(iType)-m, 'Color', 'k');
end
legend({'Unregularized LR', 'Regularized LR', 'Exponential', 'Linear'}, 'location', 'best');
ylabel('Relative log-likelihood');
title('Noise Condition');
set(gca, 'XTick', 1:N+1, 'XTickLabel', [shortnames, {'Combined'}]);
xtickangle(45);

subplot(2,1,2);
hold on;
p1 = bar([meanDiffLLRatio'; meanDiffLLRatioAll']);
drawnow;
for iSubject=1:N
    for iType=1:4
        c = p1(iType).FaceColor;
        x = p1(iType).XData(iSubject) + p1(iType).XOffset;
        m = meanDiffLLRatio(iType, iSubject);
        errorbar(x, m, m-loLLRatio(iType, iSubject), hiLLRatio(iType, iSubject)-m, 'Color', 'k');
    end
end
for iType=1:4
    c = p1(iType).FaceColor;
    x = N + 1 + p1(iType).XOffset;
    m = meanDiffLLRatioAll(iType);
    errorbar(x, m, m-loLLRatioAll(iType), hiLLRatioAll(iType)-m, 'Color', 'k');
end
ylabel('Relative log-likelihood');
title('Ratio Condition');
set(gca, 'XTick', 1:N+1, 'XTickLabel', [shortnames, {'Combined'}]);
xtickangle(45);

%% Figure 4 - insufficiency of ITB model

fig4 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

beta_range = [-0.5 .15];
pk_hprs = [1 0 0];

% Case 1: no bound, no leak (aka ideal observer)
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'bound', inf, 'gamma', 0);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 1, parula);
figureToPanel(pk_fig, fig4, 5, 3, 2);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 3, cmap);

% Case 2: bound but no leak (Kiani et al 2008)
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'bound', .8, 'gamma', 0);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 4, parula);
figureToPanel(pk_fig, fig4, 5, 3, 5);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 6, cmap);

% Case 3: leak but no bound
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'bound', inf, 'gamma', .1);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 7, parula);
figureToPanel(pk_fig, fig4, 5, 3, 8);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 9, cmap);

% Case 4: both leak and bound
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'bound', .8, 'gamma', .1);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 10, parula);
figureToPanel(pk_fig, fig4, 5, 3, 11);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 12, cmap);

% Case 5: like previous but now gamma changes with category info
params = Model.newModelParams('model', 'itb', 'trials', 10000, 'bound', .8, 'gamma_min', 0, 'gamma_max', .5);
[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, pk_hprs, 'beta', beta_range, flipud(sens_cat_pts));
figureToPanel(cs_fig, fig4, 5, 3, 13, parula);
figureToPanel(pk_fig, fig4, 5, 3, 14);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD);
figureToPanel(fig, fig4, 5, 3, 15, cmap);

%% Figure 5 - model results

fig5 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% Ideal observer performance
params = Model.newModelParams('model', 'ideal', 'trials', 10000);
disp('Loading/Running ideal observer');
[~, ideal_fig] = Model.plotCategorySensorySpace(ps, ps, params);
figureToPanel(ideal_fig, fig5, 2, 4, 5, parula);

% Sampling model with gamma = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'samples', 5);
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
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
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
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', hi_gamma, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'samples', 5);

% >> Uncomment for variational model <<
% lo_gamma = 0.1;
% hi_gamma = 0.5;
% params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', hi_gamma, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'step_size', 0.05);

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
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-0.4 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 1, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 2);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 3, cmap);

% Same with gamma = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.43 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 4, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 5);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 6, cmap);

% Same with gamma = 0.2
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.2, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.4 .23]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 7, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 8);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 9, cmap);

% --- VB-CZX ---

% Replicate vb model results above (gamma = 0)
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-0.5 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 10, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 11);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 12, cmap);

% Same with gamma = 0.1
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 13, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 14);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 15, cmap);

% Same with gamma = 0.2
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.2, 'temperature', 0.1, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.26 .2]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 16, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 17);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 18, cmap);

%% Model fits

init_model_params = Model.newModelParams('model', 'is', ...
    'gamma', 0.1, ...
    'samples', 5, ...
    'updates', 2, ...
    'lapse', 0.01);

nInner = 10;

fields = {'prior_C', 'gamma', 'samples', 'lapse', 'sensor_noise'};
LB = [0 0 1 0 0];
UB = [1 1 100 1 5];
PLB = LB;
PUB = UB;

x0 = [.5 .1 5 .01 1];

vbmc_options = vbmc('defaults');
vbmc_options.UncertaintyHandling = 'yes';

% for iSubject=1:length(noiseSubjects)
%     memo_name = ['VBMC-noise-' noiseSubjects{iSubject} '.mat'];
%     
%     SubjectData = LoadAllSubjectData(noiseSubjects{iSubject}, NOISE_PHASE, DATADIR);
%     ThresholdData = GaborThresholdTrials(SubjectData, NOISE_PHASE, 0.17, 0);
%     likefn_args = {fields, init_model_params, ThresholdData, nInner};
%     [VP, ~, ~, extflag] = LoadOrRun(@vbmc, ...
%         [{@Fitting.subjectDataLogLikelihood, x0, LB, UB, PLB, PUB, vbmc_options} likefn_args], ...
%         fullfile(MEMODIR, memo_name));
%     
%     savedir = fullfile('demo-fit', noiseSubjects{iSubject});
%     if ~exist(savedir, 'dir'), mkdir(savedir); end
%     
%     Xsamp = vbmc_rnd(VP, 1e5);
%     [fig, ax] = cornerplot(Xsamp, fields, [], [LB; UB]);
%     statuses = {'not converged', 'converged'};
%     suptitle(sprintf('%s :: %s', noiseSubjects{iSubject}, statuses{extflag+1}));
%     saveas(fig, fullfile(savedir, 'cornerplot-noise.fig'));
%     close(fig);
%     
%     fig = visModelFitPsychometric(ThresholdData, NOISE_PHASE, VP, fields, init_model_params, nInner, 10);
%     saveas(fig, fullfile(savedir, 'psychometric-noise-subthreshold.fig'));
%     close(fig);
%     
%     fig = visModelFitPsychometric(SubjectData, NOISE_PHASE, VP, fields, init_model_params, nInner, 10);
%     saveas(fig, fullfile(savedir, 'psychometric-noise.fig'));
%     close(fig);
% end

for iSubject=1:length(ratioSubjects)
    memo_name = ['VBMC-ratio-' ratioSubjects{iSubject} '.mat'];
    
    SubjectData = LoadAllSubjectData(ratioSubjects{iSubject}, RATIO_PHASE, DATADIR);
    ThresholdData = GaborThresholdTrials(SubjectData, RATIO_PHASE, 0.6, 0.4);
    likefn_args = {fields, init_model_params, ThresholdData, nInner};
    [VP, ~, ~, extflag] = LoadOrRun(@vbmc, ...
        [{@Fitting.subjectDataLogLikelihood, x0, LB, UB, PLB, PUB, vbmc_options} likefn_args], ...
        fullfile(MEMODIR, memo_name));
    
    savedir = fullfile('demo-fit', ratioSubjects{iSubject});
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    Xsamp = vbmc_rnd(VP, 1e5);
    [fig, ax] = cornerplot(Xsamp, fields, [], [LB; UB]);
    statuses = {'not converged', 'converged'};
    suptitle(sprintf('%s :: %s', ratioSubjects{iSubject}, statuses{extflag+1}));
    saveas(fig, fullfile(savedir, 'cornerplot-ratio.fig'));
    close(fig);
    
    fig = visModelFitPsychometric(ThresholdData, RATIO_PHASE, VP, fields, init_model_params, nInner, 10);
    saveas(fig, fullfile(savedir, 'psychometric-ratio-subthreshold.fig'));
    close(fig);
    
    fig = visModelFitPsychometric(SubjectData, RATIO_PHASE, VP, fields, init_model_params, nInner, 10);
    saveas(fig, fullfile(savedir, 'psychometric-ratio.fig'));
    close(fig);
end

%% Helper function for figure layout

function ax_copy = figureToPanel(figSource, figDest, subM, subN, subI, cmap)
margin = 0.02;
widths = (1 - (subN + 1) * margin) / subN;
heights = (1 - (subM + 1) * margin) / subM;
[col, row] = ind2sub([subN, subM], subI);
left = margin + (col-1) * (widths + margin);
top = margin + (row-1) * (heights + margin);
bottom = 1 - top - heights;
figure(figSource);
ax = gca;
ax_copy = copyobj(ax, figDest);
ax_copy.Position = [left bottom widths heights];
ax_copy.Parent = figDest;
close(figSource);
if exist('cmap', 'var')
    colormap(ax_copy, cmap);
    colorbar('peer', ax_copy);
end
drawnow;
end