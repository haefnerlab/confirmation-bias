function [fig, model_names] = ModelComparisonBySubject(subjectIds, phases, group_by, datadir)
%% Setup
if nargin < 3, group_by = 'subject'; end % or 'model'
if nargin < 4, datadir = fullfile(pwd, '..', 'PublishData'); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

if ~iscell(subjectIds), subjectIds = {subjectIds}; end
if ~iscell(phases), phases = {phases}; end
% Names copied from @ModelComparison
model_names = {'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split', 'ideal'};

nSubjects = length(subjectIds);
nPhases = length(phases);
nModels = length(model_names);

%% Run model fitting/scoring

aic = nan(nSubjects, nModels, nPhases);
aic_err = nan(nSubjects, nModels, nPhases);
parfor ii=1:numel(aic)
    [iSubject, iModel, iPhase] = ind2sub([nSubjects nModels nPhases], ii);

    kernel_kappa = 0.16;
    [base_params, sigs, choices] = GetSubjectDataForFitting(subjectIds{iSubject}, kernel_kappa, [], datadir);

    % Phase may be 1 for ratio, 2 for noise, or [1 2] for both
    if length(phases{iPhase}) == 1
        this_params = base_params(iPhase);
        sigs = sigs{phases{iPhase}};
        choices = choices{phases{iPhase}};

        % If single-phase, skip all 'split' models
        if contains(model_names{iModel}, 'split'), continue; end
    else
        this_params = base_params;
    end

    prefix = [subjectIds{iSubject} '-' num2str(kernel_kappa) '-' num2str(phases{iPhase})];
    [aic(ii), aic_err(ii)] = ModelComparison(this_params, sigs, choices, true, prefix, model_names{iModel});
end

%% Plot result
for iPhase=1:nPhases
    fig = figure; hold on;
    if startsWith(group_by, 's')
        h = bar(aic(:,:,iPhase)); drawnow;
        for s=1:nSubjects
            errorbar(h(s).XData+h(s).XOffset, aic(s,:,iPhase), aic_err(s,:,iPhase), 'ok');
        end
        legend(model_names, 'Location', 'best');
        set(gca, 'XTick', 1:nSubjects, 'XTickLabel', subjectIds);
    elseif startsWith(group_by, 'm')
        h = bar(aic(:,:,iPhase)'); drawnow;
        for m=1:nModels
            errorbar(h(m).XData+h(m).XOffset, aic(:,m,iPhase), aic_err(:,m,iPhase), 'ok');
        end
        legend(subjectIds, 'Location', 'best');
        set(gca, 'XTick', 1:nModels, 'XTickLabel', model_names);
    end
    grid on;
    ylabel('AIC');
    title(['Phase: ' num2str(phases{iPhase})]);
end
end