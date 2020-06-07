function [samples, fields, loglike, diagnostics, params, sigs, choices] = GetITBPosteriorSamples(subjectOrModel, phaseStr, nSamples, datadir, memodir)
%GETITBPOSTERIORSAMPLES helper function that wraps fitting posterior over the 'generic' ITB model.

if nargin < 4, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 5, memodir = fullfile(datadir, '..', 'Precomputed'); end

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
else
    kernel_kappa = 0.16;
    [params, sigs, choices] = GetSubjectDataForFitting(subjectOrModel, kernel_kappa, [], datadir);
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

distribs = Fitting.defaultDistributions(fields, false);

% Set params.model to 'itb-int' to use deterministic likelihoods from numerical integration over
% internal noise.
for iPhase=1:length(params)
    params(iPhase).model = 'itb-int';
end

% Copied from +Fitting.fitModelBADS and +Fitting.fitModelMH : construct a UID from hash of model,
% fields, input signals, and choices. Agnostic to whether this is simulated or subject data.
if iscell(sigs)
    allsigs = vertcat(sigs{:});
    allchoices = vertcat(choices{:});
    input_id = string2hash([params(1).model, strjoin(fields), num2str([allsigs(:)' allchoices'])]);
else
    input_id = string2hash([params(1).model, strjoin(fields), num2str([sigs(:)' choices'])]);
end

if ~exist('mh-checkpoints', 'dir'), mkdir('mh-checkpoints'); end
chkpt = fullfile('mh-checkpoints', sprintf('%X', input_id));

[samples, accept, smpl_info] = Fitting.sampleModelMH(...
    sigs, choices, params, nSamples, distribs, 0, 0, 1, chkpt);

loglike = smpl_info.loglike(:);

%% Once 'nSamples' drawn, get some diagnostics
diagnostics.tbl = MCMCDiagnostics(samples', fields(:), 500);
diagnostics.accept = accept;
end