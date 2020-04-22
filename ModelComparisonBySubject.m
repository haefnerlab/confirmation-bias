function [fig, model_names] = ModelComparisonBySubject(subjectIds, datadir)
%% Setup
if nargin < 2, datadir = fullfile(pwd, '..', 'PublishData'); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

if ~iscell(subjectIds), subjectIds = {subjectIds}; end

% Note that for legacy reasons, subject data in 'both' conditions are concatenated in {hslc, lshc}
% order, which is the opposite of how we order model data
phases = {2, 1, [1 2]};
phase_names = {'LSHC', 'HSLC', 'both'};

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
for iPhase=1:nPhases
    for iSubject=1:nSubjects
        ax = subplot(nPhases, nSubjects, (iPhase-1)*nSubjects+iSubject); hold on;
        h = bar(ax, [1 nan], [aic(iSubject,:,iPhase); nan(1, nModels)]);
        drawnow;
        errorbar(ax, h(1).XData(1)+[h.XOffset], aic(iSubject,:,iPhase), 3*aic_err(iSubject,:,iPhase), '.k');
        
        if iPhase == nPhases && iSubject == nSubjects, legend(model_names, 'Location', 'best'); end
        if iSubject == 1, ylabel('AIC'); end
        grid on;
        set(ax, 'XTick', []);
        ylim([min(min(aic(iSubject,:,iPhase), [], 2)) inf]-10);
        title([short_names{iSubject} ' ' phase_names{iPhase}]);
    end
end
end