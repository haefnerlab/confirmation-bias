function [samples, accept, smpl_info] = sampleModelMH(signals, choices, base_params, n_samples, distributions, ibs_repeats, sample_burnin, sample_thin, chkpt)
if nargin < 6 || isempty(ibs_repeats), ibs_repeats = 5; end
if nargin < 7 || isempty(sample_burnin), sample_burnin = 100; end
if nargin < 8 || isempty(sample_thin), sample_thin = 10; end
if nargin < 9, chkpt = ''; end

fields = fieldnames(distributions);
for iField=1:length(fields)
    if ~isfield(distributions.(fields{iField}), 'logproppdf')
        distributions.(fields{iField}).logproppdf = ...
            @(x1, x2) log(distributions.(fields{iField}).proppdf(x1, x2));
    end
end

% If ibs_repeats is set to 0, it indicates that we should use the deterministic/exact likelihood
% from @Fitting.choiceModelLogProb (which is biased if stochastic) rather than the stochastic but
% unbiased estimator @Fitting.choiceModelLogProbIBS
if ibs_repeats == 0
    if Model.isStochastic(base_params)
        warning('IBS is recommended for stochastic models! Set ibs_repeats > 0!');
    end
    log_post = @(smpl) Fitting.choiceModelLogProb(...
        Fitting.setParamsFields(base_params, fields, smpl), distributions, signals, choices, 1);
else
    if ~Model.isStochastic(base_params)
        warning('IBS is not recommended for deterministic models! Set ibs_repeats = 0!');
    end
    log_post = @(smpl) Fitting.choiceModelLogProbIBS(...
        Fitting.setParamsFields(base_params, fields, smpl), distributions, signals, choices, [], ibs_repeats);
end

    function x = proprnd(x, fields, distribs)
        for jField=1:length(fields)
            x(jField) = distribs.(fields{jField}).proprnd(x(jField));
        end
    end

    function log_prop = logproppdf(x1, x2, fields, distribs)
        log_prop = 0;
        for jField=1:length(fields)
            log_prop = log_prop + distribs.(fields{jField}).logproppdf(x1(jField), x2(jField));
        end
    end

%% Run sampler

batch_size = 100;
n_batch = ceil(n_samples / batch_size);
batch_time = nan(1, n_batch);
net_accept = nan(1, n_batch);
samples = cell(n_batch, 1);

% Load in 'checkpointed' samples if there are any
if ~isempty(chkpt)
    [pth, ~] = fileparts(chkpt);
    chkfiles = dir([chkpt '*.mat']);
    [~, isrt] = sort({chkfiles.name});
    fprintf('MHSample load %d of %d batches from checkpoints\n', length(chkfiles), n_batch);
    for iChk=length(chkfiles):-1:1
        ld = load(fullfile(pth, chkfiles(isrt(iChk)).name));
        samples{iChk} = ld.batch_sample;
        net_accept(iChk) = ld.batch_accept;
        batch_smpl_info(iChk) = ld.batch_info;
    end
    start_batch = length(chkfiles)+1;
else
    start_batch = 1;
end

% If needed (i.e. not loading from checkpoint), draw 100 points from the prior and initialize to the
% best one
if isempty(samples{1})
    disp('MHSample init');
    for iInit=100:-1:1
        init_smpl(iInit, :) = cellfun(@(f) distributions.(f).priorrnd(1), fields);
        [init_lp(iInit), ~, init_var(iInit)] = log_post(init_smpl(iInit, :));
    end
    [~,idx] = max(init_lp+3*sqrt(init_var));
    init_smpl = init_smpl(idx, :);
end

% Run sampling
disp('MHSample run');
for iBatch=start_batch:n_batch
    if iBatch > 1
        init_smpl = samples{iBatch-1}(end, :);
        sample_burnin = 0;
    end
    
    tstart = tic;
    [batch_sample, batch_accept, batch_info] = mymhsample(init_smpl, batch_size, ...
        'logpdf', log_post, ...
        'proprnd', @(x) proprnd(x, fields, distributions), ...
        'logproppdf', @(x1, x2) logproppdf(x1, x2, fields, distributions), ...
        'burnin', sample_burnin, 'thin', sample_thin);
    batch_time(iBatch) = toc(tstart);
    
    samples{iBatch} = batch_sample;
    net_accept(iBatch) = batch_accept;
    batch_smpl_info(iBatch) = batch_info;

    t_remain = (n_batch - iBatch) * nanmean(batch_time);
    fprintf('MH Sample :: Batch %03d of %03d\tAccept=%.1f%%\tETA=%.1fs\n', iBatch, n_batch, 100*nanmean(net_accept), t_remain);

    if ~isempty(chkpt)
        savename = [chkpt sprintf('-batch%04d.mat', iBatch)];
        disp(['Saving to ' savename]);
        save(savename, 'batch_sample', 'batch_accept', 'batch_info');
    end
end

samples = vertcat(samples{:});
accept = nanmean(net_accept);
smpl_info.logpdf = vertcat(batch_smpl_info.logpdf);
smpl_info.loglike = vertcat(batch_smpl_info.loglike);
smpl_info.variance = vertcat(batch_smpl_info.variance);

end