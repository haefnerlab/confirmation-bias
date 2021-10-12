function [samples, fields, tbl, params, sigs, choices] = GetStanPosteriorSamples(subjectOrModel, phaseStr, nSamples, nChains, datadir, memodir)
%GETITBPOSTERIORSAMPLES helper function that wraps fitting posterior over the 'generic' ITB model.

if nargin < 4, nChains = 1; end
if nargin < 5, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 6, memodir = fullfile(datadir, '..', 'Precomputed'); end

switch lower(phaseStr)
case 'lshc'
    phaseNo = 2;
    phaseStr = 'LSHC';
    fields = {'prior_C', 'lapse', 'noise', 'temperature', 'bound', 'neggamma'};
case 'hslc'
    phaseNo = 1;
    phaseStr = 'HSLC';
    fields = {'prior_C', 'lapse', 'noise', 'temperature', 'bound', 'neggamma'};
case 'both'
    phaseNo = [1 2];
    phaseStr = 'both';
    % Note: untying all parameters except prior, lapse, and temperature!
    fields = {'prior_C', 'lapse', 'temperature', 'noise_1', 'noise_2', 'bound_1', 'bound_2', ...
        'neggamma_1', 'neggamma_2'};
otherwise
    error('phaseStr must be one of ''HSLC'', ''LSHC'', or ''both'' (case-insensitive).');
end

% Load model simulation OR subject data into a common format
if startsWith(subjectOrModel, 'IS', 'IgnoreCase', true) || startsWith(subjectOrModel, 'ITB', 'IgnoreCase', true)
    [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {subjectOrModel, phaseNo}, ...
        fullfile(memodir, ['gt-sim-' subjectOrModel '-' phaseStr '.mat']));
    % prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phaseStr)];
    fit_scale = false;
    
    % Ensure cell format for consistency
    if ~iscell(sigs)
        sigs = {sigs};
        choices = {choices};
    end
else
    kernel_kappa = 0.16;
    [params, sigs, choices] = GetSubjectDataForFitting(subjectOrModel, kernel_kappa, [], true, datadir);
    params = params(phaseNo);
    sigs = sigs(phaseNo);
    choices = choices(phaseNo);
    
    % prefix = [subjectOrModel '-' num2str(kernel_kappa) '-' num2str(phaseNo)];
    fit_scale = true;
end

if fit_scale
    if length(phaseNo) == 1
        fields = [fields {'signal_scale'}];
    else
        % Untie scales between LSHC and HSLC fits
        fields = [fields {'signal_scale_1', 'signal_scale_2'}];
    end
end

sample_dir = 'stan-samples';
if ~exist(sample_dir, 'dir'), mkdir(sample_dir); end

data = struct('trials', length(choices{1}), 'frames', 10, 'signals', sigs{1}, 'choices', double(choices{1} == +1));
fit = stanWrapper('file', '+Fitting/functional_integrator.stan', 'data', data, ...
    'iter', nSamples, 'warmup', 500, 'chains', nChains, 'working_dir', sample_dir);
samples = fit.extract();

% Convert back to chains
samples = mat2cell(samples, size(samples,1)/nChains*ones(nChains,1), size(samples,2));

%% Once 'nSamples' drawn, get some diagnostics
[~,tbl] = fit.print();
end