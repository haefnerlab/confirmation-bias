function plotITBPM(subjectId, fitPhase, phaseNo, datadir, memodir)
if nargin < 3, phaseNo = 1; end
if nargin < 4, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 5, memodir = fullfile(datadir, '..', 'Precomputed'); end

[samples, fields] = GetITBPosteriorSamples(subjectId, fitPhase, 0, true, 1:12, datadir, memodir);
samples = vertcat(samples{:});

% Load model simulation OR subject data into a common format
if startsWith(subjectId, 'IS', 'IgnoreCase', true) || startsWith(subjectId, 'ITB', 'IgnoreCase', true)
    [params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {subjectId, phaseNo}, ...
        fullfile(memodir, ['gt-sim-' subjectId '-' phaseStr '.mat']));
    
    % Ensure cell format for consistency
    if ~iscell(sigs)
        sigs = {sigs};
        choices = {choices};
    end
else
    kernel_kappa = 0.16;
    [params, sigs, choices] = GetSubjectDataForFitting(subjectId, kernel_kappa, [], false, datadir);
    params = params(phaseNo);
    sigs = sigs(phaseNo);
    choices = choices(phaseNo);
end
sigs = sigs{phaseNo};
choices = choices{phaseNo};

trial_sigs = mean(sigs, 2);
extreme_sig = max(abs(trial_sigs));
sig_bin_edges = linspace(-extreme_sig, +extreme_sig, 13);
sig_bin_ctrs = (sig_bin_edges(1:end-1) + sig_bin_edges(2:end))/2;

bins = discretize(trial_sigs, sig_bin_edges);
for ibin=length(sig_bin_ctrs):-1:1
    [subjectPData(ibin,:), subjectEData(ibin,:)] = binofit(sum(choices(bins == ibin) == +1), sum(bins == ibin));
end

smpl_thin = 1:50:size(samples,1);
for isamp=1:length(smpl_thin)
    this_params = Fitting.setParamsFields(params, fields, samples(smpl_thin(isamp),:));
    this_params = this_params(phaseNo);
    this_params.model = 'itb';
    results = Model.runVectorized(this_params, sigs/this_params.signal_scale);
    for ibin=length(sig_bin_ctrs):-1:1
        modelFracChoice(isamp,ibin) = mean(results.choices(bins == ibin) == +1);
    end
end

[modelP, modelLo, modelHi] = meanci(modelFracChoice, 0.95);

figure; hold on;
[hl, hp] = boundedline(sig_bin_ctrs', modelP', [modelP-modelLo; modelHi-modelP]', 'alpha');
he = errorbar(sig_bin_ctrs, subjectPData, subjectPData-subjectEData(:,1), subjectEData(:,2)-subjectPData, '.k');
xlabel('mean signal');
ylabel('% chose left');
grid on;
legend([hl, hp, he], {'Model mean', 'Model 95% CI', 'data'});

end