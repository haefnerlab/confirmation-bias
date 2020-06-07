function [aic, aic_err, ll, ll_err, ntrials, fits, phase_names, model_info] = ModelComparisonBySubject(subjectIds, datadir, memodir)
if nargin < 2, per_trial = false; end
if nargin < 3, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 4, memodir = fullfile(datadir, '..', 'Precomputed'); end
if ~exist(memodir, 'dir'), mkdir(memodir); end

%% Setup

if ~iscell(subjectIds), subjectIds = {subjectIds}; end

% Note that for legacy reasons, subject data in 'both' conditions are concatenated in {hslc, lshc}
% order, which is the opposite of how we order model data
phases = {1, 2, [1 2]};
phase_names = {'HSLC', 'LSHC', 'both'};

% All names from @ModelComparison
model_names = {'ideal', 'is', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split'};

nSubjects = length(subjectIds);
nPhases = length(phases);
nModels = length(model_names);

%% Load all fits
aic = nan(nSubjects, nModels, nPhases);
parfor ii=1:numel(aic)
ll = nan(nSubjects, nModels, nPhases);
ll_err = nan(nSubjects, nModels, nPhases);
ntrials = zeros(nSubjects, nModels, nPhases);
model_info = repmat(struct('name', '', 'fields', {{}}, 'color', [0 0 0], 'type', '', 'plotname', ''), size(aic));
fits = cell(nSubjects, nModels, nPhases);
    [iSubject, iModel, iPhase] = ind2sub([nSubjects nModels nPhases], ii);

    kernel_kappa = 0.16;
    [base_params, sigs, choices] = GetSubjectDataForFitting(subjectIds{iSubject}, kernel_kappa, [], datadir);

    % Phase may be 1 for ratio, 2 for noise, or [1 2] for both
    if length(phases{iPhase}) == 1
        this_params = base_params(iPhase);
        sigs = sigs{phases{iPhase}};
        choices = choices{phases{iPhase}};

        ntrials(ii) = length(choices);

        % If single-phase, skip all 'split' models
        if contains(model_names{iModel}, 'split'), continue; end
    else
        this_params = base_params;
        ntrials(ii) = sum(cellfun(@length, choices));
    end

    prefix = [subjectIds{iSubject} '-' num2str(kernel_kappa) '-' num2str(phases{iPhase})];
    [aic(ii), ll(ii), ll_err(ii), model_info(ii), fits{ii}] = ...
        ModelComparison(this_params, sigs, choices, true, prefix, model_names{iModel}, true, memodir);
end

% since aic=-2*LL+adjustment, error in AIC is 2*err in LL
aic_err = 2*ll_err;

end