function [fig_slopes, fig_explain] = BetaExplainedBySubject(subjectId, datadir, memodir)
if nargin < 2, datadir = fullfile('..', 'PublishData'); end
if nargin < 3, memodir = fullfile(datadir, '..', 'Precomputed'); end

kernel_kappa = 0.16;

fit_phases = {2, 1, [1 2], [1 2]};
test_phases = {2, 1, 2, 1};
nPhases = length(fit_phases);
phase_names = {'LSHC', 'HSLC', 'both', 'both'};
phase_titles = {'LSHC', 'HSLC', 'Both - LSHC', 'Both - HSLC'};
stair_var = {'noise', 'true_ratio', 'noise', 'true_ratio'};

% Names copied from @ModelComparison
model_names = {'itb', 'itb-gamma', 'is'};
model_all_fields = {
    {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'neggamma', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'gamma', 'samples'}};
model_test_fields = {{'bound','gamma','noise'}, {'bound','gamma','noise'}, {'samples','gamma'}};
nModels = length(model_names);

%% Load fits

fig_slopes = figure;
for iPhase=nPhases:-1:1
    SubjectData = LoadAllSubjectData(subjectId, test_phases{iPhase}, datadir);
    [~, sigs, choices] = GetSubjectDataForFitting(subjectId, kernel_kappa, [], datadir);
    sigs = sigs{test_phases{iPhase}};
    choices = choices{test_phases{iPhase}};
    
    if test_phases{iPhase} == 2
        [floor, thresh] = GaborAnalysis.getThresholdWindow(SubjectData, test_phases{iPhase}, 0.5, 0.7, memodir);
    else
        floor = 0.4;
        thresh = 0.6;
    end
    memo_name = ['Boot-ExpPK-' num2str(kernel_kappa) '-' stair_var{iPhase} '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
    [~, ~, ~, ~, ~, abb{iPhase}] = LoadOrRun(@BootstrapExponentialWeightsGabor, {sigs, choices, 500, true}, ...
        fullfile(memodir, memo_name));
    
    fit_prefix = [subjectId '-' num2str(kernel_kappa) '-' num2str(fit_phases{iPhase})];
    for iModel=nModels:-1:1
        mname = model_names{iModel};
        if length(fit_phases{iPhase})==2 && contains(mname, 'itb')
            mname = [mname '-split'];
        end
        
        [thefit, complete] = GetFit(subjectId, phase_names{iPhase}, mname, true, datadir, memodir);
        if length(thefit) == 2
            thefit = thefit(test_phases{iPhase});
        end
        if ~isempty(thefit)
            if complete
                % If complete, we can run and store the final pk-beta-explained result
                memo_file_beta = fullfile(memodir, ['pk-beta-explained-' fit_prefix '-' num2str(test_phases{iPhase}) '-' mname '.mat']);
                [betaExplained{iPhase, iModel}, betaErr{iPhase, iModel}, ablation{iPhase, iModel}] = LoadOrRun(@PKShapeExplained, ...
                    {sigs, thefit, model_test_fields{iModel}, model_all_fields{iModel}}, memo_file_beta);
            else
                % If loaded an interim evaluation before fitting complete, we don't want to cache
                % the result
                [betaExplained{iPhase, iModel}, betaErr{iPhase,iModel}, ablation{iPhase, iModel}] = PKShapeExplained(...
                    sigs, thefit, model_test_fields{iModel}, model_all_fields{iModel});
            end
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
            
            bar(1:length(betaExplained{iPhase,iModel}), betaExplained{iPhase,iModel});
            errorbar(1:length(betaExplained{iPhase,iModel}), betaExplained{iPhase,iModel}, betaErr{iPhase,iModel}, 'ok');
            
            [meanTru, loTru, hiTru] = meanci(abb{iPhase}(:,2));
            h = bar(length(betaExplained{iPhase,iModel})+1, meanTru);
            h.FaceColor = [1 0 0];
            errorbar(length(betaExplained{iPhase,iModel})+1, meanTru, meanTru-loTru, hiTru-meanTru, 'ok');
            
            abl = ablation{iPhase, iModel};
            allfields = unique(horzcat(abl{:}));
            names = cellfun(@(a) strjoin(setdiff(allfields, a), '+'), abl, 'uniformoutput', false);
            names{cellfun(@isempty, names)} = 'null';
            names = [names {'true'}];
            set(gca, 'XTick', 1:length(names), 'XTickLabel', names);
            xtickangle(60);
            
            mname = model_names{iModel};
            if length(fit_phases{iPhase})==2 && contains(mname, 'itb')
                mname = [mname '-split'];
            end
            title([phase_titles{iPhase} ' [' mname ']']);
            
            ylim([-.4 .4]);
        end
    end
end
sgtitle(subjectId);