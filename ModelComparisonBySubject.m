function [fig, model_names] = ModelComparisonBySubject(subjectIds, phases, group_by, datadir)
%% Setup
if nargin < 3, group_by = 'subject'; end % or 'model'
if nargin < 4, datadir = fullfile(pwd, '..', 'PublishData'); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

if ~iscell(subjectIds), subjectIds = {subjectIds}; end
if ~iscell(phases), phases = {phases}; end
% Names copied from @ModelComparison
model_names = {'ideal', 'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split'};

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
short_names = cellfun(@(s) regexprep(s, '[\w-]+(subject\d+)', '$1'), subjectIds, 'uniformoutput', false);
fig = figure;
baseline = aic(:,strcmpi(model_names, 'ideal'),:);
for iPhase=1:nPhases
    ax = subplot(nPhases, 1, iPhase);
    hold on;
    if startsWith(group_by, 's')
        valid = ~all(isnan(aic(:,:,iPhase)), 1);
        thisaic = aic(:,valid,:);
        h = bar(ax, thisaic(:,:,iPhase)-baseline(:,1,iPhase)); drawnow;
        for m=1:length(h)
            errorbar(ax, h(m).XData+h(m).XOffset, thisaic(:,m,iPhase)-baseline(:,1,iPhase), aic_err(:,m,iPhase), 'o', 'Color', h(m).FaceColor/3);
        end
        legend(model_names(valid), 'Location', 'best');
        set(ax, 'XTick', 1:nSubjects, 'XTickLabel', short_names);
    elseif startsWith(group_by, 'm')
        h = bar(ax, aic(:,:,iPhase)'-baseline(:,:,iPhase)'); drawnow;
        for s=1:length(h)
            errorbar(ax, h(s).XData+h(s).XOffset, aic(s,:,iPhase)-baseline(:,:,iPhase), aic_err(s,:,iPhase), 'o', 'Color', h(s).FaceColor/3);
        end
        legend(subjectIds, 'Location', 'best');
        set(ax, 'XTick', 1:nModels, 'XTickLabel', model_names);
    end
    grid on;
    ylabel(ax, 'relative AIC - ideal');
    title(['Phase: ' num2str(phases{iPhase})]);
end
end