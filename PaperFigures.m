%% Setup

clear;
RATIO_PHASE = 1;
NOISE_PHASE = 2;
THRESHOLD = 0.7;
DATADIR = fullfile('..', 'RawData');

%% Figure 1

% Nothing to do - conceptual figures only created in vector graphics editor

%% Figure 2 - data analysis

ratioSubjects = arrayfun(@(i) sprintf('bpgFinaltest-subject%02d', i), setdiff(1:15, [1 5 12]), 'UniformOutput', false);
noiseSubjects = arrayfun(@(i) sprintf('bpgFinaltest-subject%02d', i), setdiff(1:15, [1 5]), 'UniformOutput', false);

% informed is the opposite of naive
informedSubjects = arrayfun(@(i) sprintf('bpgFinaltest-subject%02d', i), [7 8 9], 'UniformOutput', false);

naiveRatioSubjects = setdiff(ratioSubjects, informedSubjects);
naiveNoiseSubjects = setdiff(noiseSubjects, informedSubjects);
naiveBothSubjects = intersect(naiveRatioSubjects, naiveNoiseSubjects);

% Re-order so that informed subjects are at the end
ratioSubjects = [naiveRatioSubjects informedSubjects];
noiseSubjects = [naiveNoiseSubjects informedSubjects];
bothSubjects = [intersect(naiveRatioSubjects, naiveNoiseSubjects) informedSubjects];

% Compute population-level psycho metrics
N = length(bothSubjects);
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
fprintf('Noise condition mean threshold = %f +/- %f\n', mean(threshold), sem(threshold));
fprintf('Noise condition mean %%correct at zero signal = %f +/- %f\n', mean(100*pc_zero_signal), sem(100*pc_zero_signal));
fprintf('Ratio condition mean %%correct at 6:4 = %f +/- %f\n', mean(100*pc_60_40_adjusted), sem(100*pc_60_40_adjusted));

GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [1 2], 'exponential');
GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [1 2], 'linear');
GaborAnalysis.DeltaPK(bothSubjects, [1 2], false);

%% Figure 3 - sampling model results

fig3 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% Ideal observer performance
params = Model.newModelParams('model', 'ideal', 'trials', 10000);
[~, ideal_fig] = Model.plotCategorySensorySpace(ps, ps, params);
figureToPanel(ideal_fig, fig3, 2, 4, 5, parula);

% Model performance with gamma = 0.
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-0.5 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [.99 .59; .83 .65; .69 .79; .67 .91];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig3, 2, 4, 2, parula);
figureToPanel(pk_fig, fig3, 2, 4, 3);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, fig3, 2, 4, 4, cmap);

% Same plots with gamma = 0.1 showing emergence of recency effects
params.gamma = 0.1;
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.61; 0.85 0.65; 0.73 0.73; 0.67 0.83; 0.65 0.95];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig3, 2, 4, 6, parula);
figureToPanel(pk_fig, fig3, 2, 4, 7);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, fig3, 2, 4, 8, cmap);

%% Supplemental figure model results (noise term + VB model)

figModelSupp = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% --- SAMPLING ---

% Replicate sampling model results above (gamma = 0)
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-0.5 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [.99 .59; .83 .65; .69 .79; .67 .91];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 1, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 2);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 3, cmap);

% Same with gamma = 0.1, noise = 0
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.61; 0.85 0.65; 0.73 0.73; 0.67 0.83; 0.65 0.95];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 4, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 5);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 6, cmap);

% Same with gamma = 0, noise = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'noise', 0.1, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 .59; 0.81 .67; 0.69 .81; 0.65 .99];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 7, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 8);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 9, cmap);

% --- VB-CZX ---

% Replicate sampling model results above (gamma = 0)
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0, 'noise', 0, 'trials', 10000, 'updates', 5, 'step_size', 0.01);
beta_range = [-0.5 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.59; 0.81 0.65; 0.73 0.73; 0.65 0.85; 0.65 0.97];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 10, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 11);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 12, cmap);

% Same with gamma = 0.1, noise = 0
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'step_size', 0.01);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.59; 0.85 0.63; 0.73 0.71; 0.65 0.83; 0.63 0.97];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range);
figureToPanel(cs_fig, figModelSupp, 6, 3, 13, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 14);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 15, cmap);

% Same with gamma = 0, noise = 0.05
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0, 'noise', 0.05, 'trials', 10000, 'updates', 5, 'step_size', 0.01);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.69; 0.91 0.75; 0.85 0.83; 0.81 0.93; 0.79 0.99];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range);
figureToPanel(cs_fig, figModelSupp, 6, 3, 16, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 17);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 18, cmap);

%% Helper function for figure layout

function figureToPanel(figSource, figDest, subM, subN, subI, cmap)
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
end