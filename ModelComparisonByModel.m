function fig = ModelComparisonByModel(memodir)
gt_names = {'IS', 'ITB'};
nTruth = length(gt_names);

% Note that for legacy reasons, model data are concatenated in {lshc, hslc} order, which is the
% opposite of how we usually concatenate human data
phases = {2, 1, [2 1]};
nPhases = length(phases);
phase_names = {'LSHC', 'HSLC', 'both'};

% Names copied from @ModelComparison
model_names = {'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split', 'ideal'};
nModels = length(model_names);
%%

aic = nan(nTruth, nModels, nPhases);
aic_err = nan(nTruth, nModels, nPhases);
parfor ii=1:numel(aic)
    [iTruth, iModel, iPhase] = ind2sub([nTruth, nModels, nPhases], ii);

    [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {gt_names{iTruth}, phases{iPhase}}, ...
        fullfile(memodir, ['gt-sim-' gt_names{iTruth} '-' phase_names{iPhase} '.mat']));
    
    if contains(model_names{iModel}, 'split') && length(phases{iPhase}) < 2
        continue;
    end
    
    prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phase_names{iPhase})];
    [aic(ii), aic_err(ii)] = ModelComparison(params, sigs, choices, false, prefix, model_names(iModel));
end

%% Plot result
fig = figure;
for iPhase=1:nPhases
    ax = subplot(nPhases, 1, iPhase); hold on;
    h = bar(aic(:,:,iPhase)); drawnow;
    for s=1:nTruth
        errorbar(h(1).XData(s)+[h.XOffset], aic(s,:,iPhase), aic_err(s,:,iPhase), '.k');
    end
    if iPhase == 1, legend(model_names, 'Location', 'best'); end
    set(gca, 'XTick', 1:nTruth, 'XTickLabel', gt_names);
    grid on;
    ylim([min(min(aic(:,:,iPhase), [], 2)) inf]-50);
    ylabel('AIC');
    title(phase_names{iPhase});
end
end