%% Figure 1

% Nothing to do - conceptual figures only created in vector graphics editor

%% Figure 2 - data analysis

RATIO_PHASE = 1;
NOISE_PHASE = 2;
DATADIR = fullfile('..', 'RawData');

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
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, NOISE_PHASE, 0.5, 0.75, DATADIR);
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

% Model performance with gamma = 0. (Note that this requires user input to get points along
% threshold curve)
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-0.5 -eps]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [.99 .59; .83 .65; .69 .79; .67 .91];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig3, 2, 4, 2, parula);
figureToPanel(pk_fig, fig3, 2, 4, 3);
[~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range);
figureToPanel(fig, fig3, 2, 4, 4, cmap);

% Same plots with gamma = 0.1 showing emergence of recency effects
params.gamma = 0.1;
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts =[0.99 0.61; 0.85 0.65; 0.73 0.73; 0.67 0.83; 0.65 0.95];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig3, 2, 4, 6, parula);
figureToPanel(pk_fig, fig3, 2, 4, 7);
[~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range);
figureToPanel(fig, fig3, 2, 4, 8, cmap);

%% Figure 3 (variational bayes version)

% fig3vb = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% Ideal observer performance
params = Model.newModelParams('model', 'ideal', 'trials', 10000);
Model.plotCategorySensorySpace(ps, ps, params);

% Model performance with gamma = 0. (Note that this requires user input to get points along
% threshold curve)
params = Model.newModelParams('model', 'vb', 'gamma', 0, 'trials', 10000, 'updates', 5, 'noise', .5);
beta_range = [-5.5 0]; % min and max beta expected (to get maximum use of colorbar range)
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range);
Model.plotCSSlopes(ps, ps, params, beta_range);

% Same plots with gamma = 0.2
params.gamma = 0.2;
beta_range = [-.2 .2]; % min and max beta expected (to get maximum use of colorbar range)
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range);
Model.plotCSSlopes(ps, ps, params, beta_range);

function figureToPanel(figSource, figDest, subM, subN, subI, cmap)
figure(figSource);
ax = gca;
ax_copy = copyobj(ax, figDest);
ax_tmp = subplot(subM, subN, subI);
ax_copy.Position = ax_tmp.Position;
ax_copy.Parent = figDest;
close(figSource);
if exist('cmap', 'var')
    colormap(ax_copy, cmap);
    colorbar('peer', ax_copy);
end
end