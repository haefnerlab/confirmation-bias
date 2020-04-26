function [fig_slopes, fig_explain] = BetaExplainedBySubject(subjectId, datadir)
if nargin < 2, datadir = fullfile('..', 'PublishData'); end
memodir = fullfile(datadir, '..', 'Precomputed');

kernel_kappa = 0.16;

phases = {2, 1};
nPhases = length(phases);
phase_names = {'LSHC', 'HSLC'};
stair_var = {'noise', 'true_ratio'};

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
for iPhase=nPhases:-1:1
    SubjectData = LoadAllSubjectData(subjectId, phases{iPhase}, datadir);
    [~, sigs, choices] = GetSubjectDataForFitting(subjectId, kernel_kappa, [], datadir);
    sigs = sigs{phases{iPhase}};
    choices = choices{phases{iPhase}};
    
    prefix = [subjectId '-' num2str(kernel_kappa) '-' num2str(phases{iPhase})];
    
    if phases{iPhase} == 2
        [floor, thresh] = GaborAnalysis.getThresholdWindow(SubjectData, phases{iPhase}, 0.5, 0.7, memodir);
    else
        floor = 0.4;
        thresh = 0.6;
    end
    memo_name = ['Boot-ExpPK-' num2str(kernel_kappa) '-' stair_var{iPhase} '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
    [~, ~, ~, ~, ~, abb{iPhase}] = LoadOrRun(@BootstrapExponentialWeightsGabor, {sigs, choices, 500, true}, ...
        fullfile(memodir, memo_name));
    
    for iModel=nModels:-1:1
        memo_file_fit = fullfile(memodir, ['qrgfit-' prefix '-' model_names{iModel} '.mat']);
        if exist(memo_file_fit, 'file')
            memo_file_beta = fullfile(memodir, ['pk-beta-explained-' prefix '-' model_names{iModel} '.mat']);
            ld = load(memo_file_fit);
            
            fits = ld.results{1};
            fit_names = fieldnames(fits);
            fit_lls = cellfun(@(name) fits.(name).ll, fit_names);
            [~, iBest] = max(fit_lls);
            bestfit_params = fits.(fit_names{iBest});
            
            [betaExplained{iPhase, iModel}, ablation{iPhase, iModel}] = LoadOrRun(@PKShapeExplained, ...
                {sigs, bestfit_params, model_test_fields{iModel}, model_all_fields{iModel}}, memo_file_beta);
        else
            disp(['No dice: ' memo_file_fit]);
            betaExplained{iPhase, iModel} = [];
        end
    end
    
    figure(fig_slopes); hold on;
    histogram(abb{iPhase}(:,2), -.35:.01:.35, 'DisplayName', phase_names{iPhase});
end
legend();
xlabel('PK Slope (\beta)');
title(subjectId);

%% Plot stuff

fig_explain = figure;
for iPhase=1:nPhases
    for iModel=1:nModels
        if ~isempty(betaExplained{iPhase, iModel})
            subplot(nModels, nPhases, (iModel-1)*nPhases + iPhase);
            hold on;
            
            [meanB, loB, hiB] = meanci(squeeze(mean(betaExplained{iPhase, iModel}, 1))');
            bar(1:length(meanB), meanB);
            errorbar(1:length(meanB), meanB, meanB-loB, hiB-meanB, 'ok');
            
            [meanTru, loTru, hiTru] = meanci(abb{iPhase}(:,2));
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
            
            title([phase_names{iPhase} ' [' model_names{iModel} ']']);
            
            ylim([-.35 .35]);
        end
    end
end
sgtitle(subjectId);