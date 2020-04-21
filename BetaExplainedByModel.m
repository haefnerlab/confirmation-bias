gt_names = {'IS', 'ITB'};
nTruth = length(gt_names);

% Note that for legacy reasons, model data are concatenated in {lshc, hslc} order, which is the
% opposite of how we usually concatenate human data
phases = {2, 1};
nPhases = length(phases);
phase_names = {'LSHC', 'HSLC'};

% Names copied from @ModelComparison
model_names = {'itb', 'itb-gamma', 'is'};
model_all_fields = {
    {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'neggamma', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'gamma', 'samples'}};
model_test_fields = {{'bound','gamma'}, {'bound','gamma'}, {'samples','gamma'}};
nModels = length(model_names);

%% Load fits

fig_slopes = figure;
for iTruth=nTruth:-1:1
    for iPhase=nPhases:-1:1
        [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {gt_names{iTruth}, phases{iPhase}}, ...
            fullfile(MEMODIR, ['gt-sim-' gt_names{iTruth} '-' phase_names{iPhase} '.mat']));
        
        prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phase_names{iPhase})];
        
        memo_file_beta = fullfile(MEMODIR, ['BootPK-Exp-' prefix '.mat']);
        [~, ~, ~, ~, ~, abb{iTruth, iPhase}] = LoadOrRun(@BootstrapExponentialWeightsGabor, ...
            {sigs, choices, 500, false}, memo_file_beta);
        
        for iModel=nModels:-1:1
            memo_file_fit = fullfile(MEMODIR, ['mhfit-' prefix '-' model_names{iModel} '.mat']);
            if exist(memo_file_fit, 'file')
                memo_file_beta = fullfile(MEMODIR, ['pk-beta-explained-' prefix '-' model_names{iModel} '.mat']);
                ld = load(memo_file_fit);
                bestfit_params = ld.results{1}.gp_mle_params;
                [betaExplained{iTruth, iPhase, iModel}, ablation{iPhase, iModel}] = LoadOrRun(@PKShapeExplained, ...
                    {sigs, bestfit_params, model_test_fields{iModel}, model_all_fields{iModel}}, memo_file_beta);
            else
                disp(['No dice: ' memo_file_fit]);
                betaExplained{iTruth, iPhase, iModel} = [];
            end
        end
        
        figure(fig_slopes); hold on;
        histogram(abb{iTruth, iPhase}(:,2), -.25:.01:.25, 'DisplayName', [gt_names{iTruth} ' ' phase_names{iPhase}]);
    end
end
legend();
xlabel('PK Slope (\beta)');

%% Plot stuff

fig_explain = figure;
for iCondition=1:nTruth*nPhases
    [iTruth, iPhase] = ind2sub([nTruth nPhases], iCondition);
    for iModel=1:nModels        
        if ~isempty(betaExplained{iTruth, iPhase, iModel})            
            subplot(nModels, nTruth*nPhases, (iModel-1)*nTruth*nPhases + iCondition);
            hold on;
            
            [meanB, loB, hiB] = meanci(squeeze(mean(betaExplained{iTruth, iPhase, iModel}, 1))');
            bar(1:length(meanB), meanB);
            errorbar(1:length(meanB), meanB, meanB-loB, hiB-meanB, 'ok');
            
            [meanTru, loTru, hiTru] = meanci(abb{iTruth, iPhase}(:,2));
            h = bar(length(meanB)+1, meanTru);
            h.FaceColor = [1 0 0];
            errorbar(length(meanB)+1, meanTru, meanTru-loTru, hiTru-meanTru, 'ok');
            
            abl = ablation{iPhase, iModel};
            allfields = unique(horzcat(abl{:}));
            names = cellfun(@(a) strjoin(setdiff(allfields, a), '+'), abl, 'uniformoutput', false);
            names{cellfun(@isempty, names)} = 'null';
            names = [names {'true'}];
            set(gca, 'XTick', 1:length(names), 'XTickLabel', names);
            xtickangle(60);
            
            title([gt_names{iTruth} ' ' phase_names{iPhase} ' [' model_names{iModel} ']']);
        end
    end
end