function fig = ModelComparisonByModel(memodir)
% gt_names = {'IS', 'ITB'};
gt_names = {'ITB'};
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
    
    prefix = ['gt-' Model.getModelStringID(params, true) '-' lower(phase_names{iPhase})];
    [aic(ii), aic_err(ii)] = ModelComparison(params, sigs, choices, false, prefix, model_names(iModel));
end

%% Plot result
for iPhase=1:nPhases
    fig = figure; hold on;
    h = bar(aic(:,:,iPhase)); drawnow;
    for s=1:nTruth
        errorbar(h(s).XData+h(s).XOffset, aic(s,:,iPhase), aic_err(s,:,iPhase), 'ok');
    end
    legend(model_names, 'Location', 'best');
    set(gca, 'XTick', 1:nTruth, 'XTickLabel', gt_names);
    grid on;
    ylabel('AIC');
    title(['Phase: ' num2str(phases{iPhase})]);
end
end