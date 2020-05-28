function [bestfit, allfits, complete] = GetFit(subjectOrModel, phaseStr, modelType, bestSoFar, datadir, memodir)
if nargin < 4, bestSoFar = true; end
if nargin < 5, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 6, memodir = fullfile(datadir, '..', 'Precomputed'); end

model_info = struct(...
    'name', {'ideal', 'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split'}, ...
    'type', {'ideal', 'is', 'vb-czx', 'itb', 'itb', 'itb', 'itb'}, ...
    'fields', {{'prior_C', 'log_lapse'}, ...
    {'prior_C', 'log_lapse', 'gamma', 'samples'}, ...
    {'prior_C', 'log_lapse', 'gamma', 'step_size', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'neggamma', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'gamma_1', 'gamma_2', 'bound', 'log_noise'}, ...
    {'prior_C', 'log_lapse', 'neggamma_1', 'neggamma_2', 'bound', 'log_noise'}});

switch lower(phaseStr)
case 'lshc'
    phaseNo = 2;
    phaseStr = 'LSHC';
case 'hslc'
    phaseNo = 1;
    phaseStr = 'HSLC';
case 'both'
    phaseNo = [1 2];
    phaseStr = 'both';
otherwise
    error('phaseStr must be one of ''HSLC'', ''LSHC'', or ''both'' (case-insensitive).');
end

if startsWith(subjectOrModel, 'IS', 'IgnoreCase', true) || startsWith(subjectOrModel, 'ITB', 'IgnoreCase', true)
    % Model convention for legacy reasons: 'both' is ordered [2 1]
    phaseNo = fliplr(phaseNo);

    [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {subjectOrModel, phaseNo}, ...
        fullfile(memodir, ['gt-sim-' subjectOrModel '-' phaseStr '.mat']));
    prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phaseStr)];
    fit_scale = false;
else
    kernel_kappa = 0.16;
    [params, sigs, choices] = GetSubjectDataForFitting(subjectOrModel, kernel_kappa, [], datadir);
    params = params(phaseNo);
    sigs = sigs(phaseNo);
    choices = choices(phaseNo);
    
    prefix = [subjectOrModel '-' num2str(kernel_kappa) '-' num2str(phaseNo)];
    fit_scale = true;
end

cacheFileName = fullfile(memodir, ['badsfit-' prefix '-' lower(modelType) '.mat']);

if exist(cacheFileName, 'file') || bestSoFar
    [~, ~, ~, ~, allfits, smpls] = ModelComparison(params, sigs, choices, fit_scale, prefix, modelType, bestSoFar, memodir);
    [~, ibest] = max(cellfun(@(f) f.ll, allfits{1}));
    bestfit = allfits{1}{ibest};
    % Hack...
    complete = ~isempty(smpls{1});
else
    bestfit = [];
    allfits = {};
    complete = false;
end