function fig = ModelComparisonByModel(memodir)
gt_names = {'IS', 'ITB'};
gt_nparams = [4, 5];
nTruth = length(gt_names);

% Note that for legacy reasons, model data are concatenated in {lshc, hslc} order, which is the
% opposite of how we usually concatenate human data
phases = {2, 1, [2 1]};
nPhases = length(phases);
phase_names = {'LSHC', 'HSLC', 'both'};

% Names copied from @ModelComparison
model_names = {'ideal', 'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split'};
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

%% Add ground-truth evaluations on themselves for reference
for iTruth=1:nTruth
    for iPhase=1:length(phases)
        [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {gt_names{iTruth}, phases{iPhase}}, ...
            fullfile(memodir, ['gt-sim-' gt_names{iTruth} '-' phase_names{iPhase} '.mat']));
        prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phase_names{iPhase})];
        
        [~, gt_ll, gt_ll_var] = LoadOrRun(@Fitting.choiceModelLogProbIBS, ...
            {params, struct(), sigs, choices, [], round(10*sqrt(sum([params.trials])))}, ...
            fullfile(memodir, ['gt-ll-ibs-' prefix '.mat']));
        
        aic(iTruth, nModels+1, iPhase) = aicbic(gt_ll, gt_nparams(iTruth));
        aic_err(iTruth, nModels+1, iPhase) = sqrt(gt_ll_var);
    end
end

nModels = nModels + 1;
model_names{end+1} = 'truth';

%% Plot result
fig = figure;
for iPhase=1:nPhases
    for iTruth=1:nTruth
        ax = subplot(nPhases, nTruth, (iPhase-1)*nTruth+iTruth); hold on;
        h = bar(ax, [1 nan], [aic(iTruth,:,iPhase); nan(1, nModels)]);
        drawnow;
        errorbar(ax, h(1).XData(1)+[h.XOffset], aic(iTruth,:,iPhase), 3*aic_err(iTruth,:,iPhase), '.k');
        
        if iPhase == nPhases && iTruth == nTruth, legend(model_names, 'Location', 'best'); end
        grid on;
        set(ax, 'XTick', []);
        ylim([min(min(aic(iTruth,:,iPhase), [], 2)) inf]-50);
        ylabel('AIC');
        title([gt_names{iTruth} ' ' phase_names{iPhase}]);
    end
end
end