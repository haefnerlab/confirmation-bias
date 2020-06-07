function fig = ModelComparisonByModel(memodir)
if nargin < 1, memodir = fullfile(pwd, '..', 'Precomputed'); end

gt_names = {'IS', 'ITB'};
gt_nparams = [4, 5];
nTruth = length(gt_names);

% Note that for legacy reasons, model data are concatenated in {lshc, hslc} order, which is the
% opposite of how we usually concatenate human data
phases = {2, 1, [2 1]};
nPhases = length(phases);
phase_names = {'LSHC', 'HSLC', 'both'};

% Names copied from @ModelComparison (sans vb)
model_names = {'ideal', 'is', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split'};
nModels = length(model_names);

%%
aic = nan(nTruth, nModels, nPhases);
parfor ii=1:numel(aic)
ll_err = nan(nTruth, nModels, nPhases);
model_info = repmat(struct('name', '', 'fields', {{}}, 'color', [0 0 0], 'type', '', 'plotname', ''), size(aic));
    [iTruth, iModel, iPhase] = ind2sub([nTruth, nModels, nPhases], ii);

    [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {gt_names{iTruth}, phases{iPhase}}, ...
        fullfile(memodir, ['gt-sim-' gt_names{iTruth} '-' phase_names{iPhase} '.mat']));
    
    if contains(model_names{iModel}, 'split') && length(phases{iPhase}) < 2
        continue;
    end
    
    prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phase_names{iPhase})];
    [aic(ii), ~, ll_err(ii), model_info(ii)] = ModelComparison(params, sigs, choices, false, prefix, model_names(iModel), true, memodir);
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
        ll_err(iTruth, nModels+1, iPhase) = sqrt(gt_ll_var);
        model_info(iTruth, nModels+1, iPhase).name = 'gt';
        model_info(iTruth, nModels+1, iPhase).plotname = 'Truth';
        model_info(iTruth, nModels+1, iPhase).color = [.9 .9 .9];
    end
end

nModels = nModels + 1;
model_names{end+1} = 'truth';

% Since AIC=-2*LL+parameteradjustment, the error in AIC is 2x error in LL
aic_err = 2*ll_err;

%% Plot results for each condition separately
fig = figure;
for iPhase=1:2
    for iTruth=1:nTruth
        % Subplot order: [model 1 hslc, model 1 lshc, model 2 hslc, model 2 lshc]
        ax = subplot(2, 2*nTruth, (iTruth-1)*2+iPhase); hold on;
        h = bar(ax, [1 nan], [aic(iTruth,:,iPhase); nan(1, nModels)]);
        for b=1:length(h)
            set(h(b), 'FaceColor', model_info(iTruth,b,iPhase).color);
        end
        drawnow;
        errorbar(ax, h(1).XData(1)+[h.XOffset], aic(iTruth,:,iPhase), aic_err(iTruth,:,iPhase), '.k');
        grid on;
        set(ax, 'XTick', []);
        ylim([min(min(aic(:,:,iPhase), [], 2))-10 inf]);
        ylabel('AIC');
        title([gt_names{iTruth} ' ' phase_names{iPhase}]);
    end
end

%% Plot results on joint conditions
    
% 3 groups of models: the "fully split" models equal to sum of AIC of fits to conditions
% separately, the "fully tied" models, and the 2 "gamma split" models.
idxModelJoint = cellfun(@(nm) contains(nm, '-split') || strcmp(nm, 'truth'), model_names);
idxModelToCombine = ~idxModelJoint;

plot_names = {model_info(1,:,3).plotname};
colors = vertcat(model_info(1,:,3).color);

for iTruth=1:nTruth
    subplot(2,nTruth,nTruth+iTruth); hold on;

    % AIC on joint data of combined per-condition fits is just the sum of AIC for each condition
    aic_per = aic(iTruth, idxModelToCombine, 1) + aic(iTruth, idxModelToCombine, 2);
    aic_per_err = sqrt(aic_err(iTruth, idxModelToCombine, 1).^2 + aic_err(iTruth, idxModelToCombine, 2).^2);
    per_names = cellfun(@(nm) [nm '-per'], plot_names(idxModelToCombine), 'uniformoutput', false);
    
    aic_tied = aic(iTruth, idxModelToCombine, 3);
    aic_tied_err = aic_err(iTruth, idxModelToCombine, 3);
    tied_names = plot_names(idxModelToCombine);
    
    % Reorder results so each model's 'per condition' fit is adjacent to its 'tied' fit
    idx_stack = reshape(1:2*sum(idxModelToCombine), sum(idxModelToCombine), 2)';
    idx_stack = idx_stack(:);
    
    aic_stack = [aic_per aic_tied];
    aic_err_stack = [aic_per_err aic_tied_err];
    color_stack = [colors(idxModelToCombine,:)*2/3; colors(idxModelToCombine,:)*2/3+.33];
    names_stack = [per_names tied_names];
    
    aic_gamma = aic(iTruth, idxModelJoint, 3);
    aic_gamma_err = aic_err(iTruth, idxModelJoint, 3);
    gamma_names = plot_names(idxModelJoint);
    
    allcolors = [color_stack(idx_stack,:); colors(idxModelJoint,:)];
    
    h = bar([1 nan], [aic_stack(idx_stack) aic_gamma; nan(1, 2*sum(idxModelToCombine)+sum(idxModelJoint))]);
    for b=1:length(h)
        set(h(b), 'FaceColor', allcolors(b,:));
    end
    drawnow; % populate h.XOffset
    errorbar(h(1).XData(1) + [h.XOffset], [aic_stack(idx_stack) aic_gamma], [aic_err_stack(idx_stack) aic_gamma_err], '.k');

    grid on;
    set(gca, 'XTick', []);
    ylim([min(aic(iTruth,:,3))-10 inf]);
    ylabel('AIC');
    title([gt_names{iTruth} ' Joint']);
end

legend([names_stack(idx_stack) gamma_names], 'location', 'northeast');
end
