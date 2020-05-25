function [optim_results, exitflag, gpinfo, metadata] = fitModelBADS(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELBADS Search for maximum likelhood model using Bayesian Adaptive Direct Search
%(BADS) library.

fields = fieldnames(distribs);
nF = length(fields);

% Get a UID for the model being fit, signals, and choices so that we can restart from checkpoints
if iscell(signals)
    allsigs = vertcat(signals{:});
    allchoices = vertcat(choices{:});
    input_id = string2hash([base_params(1).model, strjoin(fields), num2str([allsigs(:)' allchoices'])]);
    nTrials = length(allchoices);
else
    input_id = string2hash([base_params(1).model, strjoin(fields), num2str([signals(:)' choices'])]);
    nTrials = length(choices);
end

%% Set up evaluation, handling stochastic vs deterministic case
if Model.isStochastic(base_params)
    disp('fitModelBADS (stochastic case)');
    
    % BADS expects noise variance on the order of 1, which requires scaling # evaluations by
    % sqrt(#trials)
    ibs_repeats = round(sqrt(nTrials));
    ibs_bailout = 100;
    
    eval_fun = @Fitting.choiceModelLogProbIBS;
    eval_args = {[], ibs_repeats, ibs_bailout};
    final_eval_args = {[], 10*ibs_repeats, 10*ibs_bailout};
else
    disp('fitModelBADS (deterministic case)');
    eval_fun = @Fitting.choiceModelLogProb;
    eval_args = {1};
    final_eval_args = {1};
end

%% Select initial points using a quasi-random grid sampled from the prior
% See @Fitting.fitModelQRG

% To gracefully save/load from partially computed results, we initialize to NaN values, (maybe) load
% from a file, then only evaluate those rows that are still NaN, with periodic checkpointing.
chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-grid.mat', input_id));
if exist(chkpt, 'file')
    ld = load(chkpt);
    grid_points = ld.grid_points;
    logpdf = ld.logpdf;
    loglike = ld.loglike;
    variance = ld.variance;
    
    fprintf('fitModelBADS :: %s\tloaded %d of %d evaluations\n', chkpt, sum(~isnan(loglike)), length(loglike));
else
    % Start with a QRG of 5000 points on the unit hypercube
    hyperQRG = lhsdesign(5000, nF);
    
    % Use points inside hypercube to get grid on actual parameters via inverse CDF of each field's prior
    grid_points = zeros(size(hyperQRG));
    for iF=1:nF
        grid_points(:,iF) = distribs.(fields{iF}).priorinv(hyperQRG(:,iF));
    end
    
    % Initialize all evaluations to NaN
    logpdf = nan(size(hyperQRG,1), 1);
    loglike = nan(size(hyperQRG,1), 1);
    variance = nan(size(hyperQRG,1), 1);
end

tstart = tic;
% Run grid evaluation in batches of 'parfor' evaluations, stopping to save results after each batch.
needs_eval = find(isnan(loglike));
batchsize = 100;
nbatch = ceil(length(needs_eval)/batchsize);
for ibatch=1:nbatch
    batch_idxs = needs_eval(1:min(batchsize, length(needs_eval)));
    parfor ii=1:length(batch_idxs)
        iVal = batch_idxs(ii);
        [batch_logpdf(ii), batch_loglike(ii), batch_variance(ii)] = ...
            eval_fun(Fitting.setParamsFields(base_params, fields, grid_points(iVal, :)), distribs, signals, choices, eval_args{:});
    end
    logpdf(batch_idxs) = batch_logpdf;
    loglike(batch_idxs) = batch_loglike;
    variance(batch_idxs) = batch_variance;
    
    % Print diagnostics and save things
    telapse = toc(tstart);
    fprintf('%s :: batch %d of %d\tETA=%.1fs\n', chkpt, ibatch, nbatch, (nbatch-ibatch)*telapse/ibatch);
    save(chkpt, 'grid_points', 'logpdf', 'loglike', 'variance');
    
    % Remove indices from 'needs_eval' that were just completed.
    needs_eval(1:min(batchsize, length(needs_eval))) = [];
end

% Sanity-check
assert(~any(isnan(loglike)));

% Sort LL to start BADS from each of the top K points
[~, sort_ll] = sort(loglike, 'descend');

%% Run BADS multiple times

LB = cellfun(@(f) distribs.(f).lb, fields)';
UB = cellfun(@(f) distribs.(f).ub, fields)';
PLB = cellfun(@(f) distribs.(f).plb, fields)';
PUB = cellfun(@(f) distribs.(f).pub, fields)';

opts = bads('defaults');
opts.NonlinearScaling = 'off';
if Model.isStochastic(base_params)
    opts.UncertaintyHandling = 'on';
    opts.NoiseScale = sqrt(mean(variance(loglike > median(loglike))));
else
    opts.UncertaintyHandling = 'off';
    opts.NoiseScale = 0;
end

nInitialRuns = 15;
parfor iRun=1:nInitialRuns
    % Use the top 'nRuns' points from the grid for initialization
    chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-search-%d.mat', input_id, iRun));
    if exist(chkpt, 'file')
        ld = load(chkpt);
        this_max_x = ld.this_max_x;
        this_nll = ld.this_nll;
        this_exitflag = ld.this_exitflag;
        this_gpinfo = ld.this_gpinfo;
        fprintf('fitModelBADS :: loading search iteration %d\n', iRun);
    else
        fprintf('fitModelBADS :: starting search iteration %d\n', iRun);
        [this_max_x, this_nll, this_exitflag, ~, ~, this_gpinfo] = bads(...
            @(x) -eval_fun(Fitting.setParamsFields(base_params, fields, x), distribs, signals, choices, eval_args{:}), ...
            grid_points(sort_ll(iRun), :), LB, UB, PLB,PUB, [], opts);
        saveWrapper(chkpt, struct('this_max_x', this_max_x, 'this_nll', this_nll, ...
            'this_exitflag', this_exitflag, 'this_gpinfo', this_gpinfo));
    end
    max_x(iRun,:) = this_max_x;
    nll(iRun) = this_nll;
    exitflag(iRun) = this_exitflag;
    gpinfo(iRun) = this_gpinfo;
end

%% Run BADS an additional number of times until confident in global min

% Extrapolate from first nInitialRuns out to get minimum number of runs such that we could be
% confident that they contain a global minimum. Since BADS internally runs 10 evaluations of the
% final point by default, an estimate of the noise around the minima is opts.NoiseScale/sqrt(10).
estMinRuns = Fitting.bootstrapMinimaRegret(nll, opts.NoiseScale/sqrt(10), 1);
fprintf('fitModelBADS :: completed %d initial runs; estimated low-regret # runs = %d\n', nInitialRuns, estMinRuns);
if estMinRuns > nInitialRuns
    % Copy of the above parfor loop, but now initializing by selecting among the bestfit points so
    % far plus some jittering
    jitteramt = (PUB-PLB)/20;
    parfor iAddRun=1:estMinRuns-nInitialRuns
        iRun = iAddRun+nInitialRuns;
        
        chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-search-%d.mat', input_id, iRun));
        if exist(chkpt, 'file')
            ld = load(chkpt);
            this_max_x = ld.this_max_x;
            this_nll = ld.this_nll;
            this_exitflag = ld.this_exitflag;
            this_gpinfo = ld.this_gpinfo;
            fprintf('fitModelBADS :: loading additional search iteration %d\n', iRun);
        else
            fprintf('fitModelBADS :: starting additional search iteration %d\n', iRun);
            start_pt = max_x(randi(nInitialRuns), :) + randn(1,nF).*jitteramt;
            start_pt = max(PLB, min(start_pt, PUB));
            [this_max_x, this_nll, this_exitflag, ~, ~, this_gpinfo] = bads(...
                @(x) -eval_fun(Fitting.setParamsFields(base_params, fields, x), distribs, signals, choices, eval_args{:}), ...
                start_pt, LB, UB, PLB, PUB, [], opts);
            saveWrapper(chkpt, struct('this_max_x', this_max_x, 'this_nll', this_nll, ...
                'this_exitflag', this_exitflag, 'this_gpinfo', this_gpinfo));
        end
        add_max_x(iAddRun,:) = this_max_x;
        add_nll(iAddRun) = this_nll;
        add_exitflag(iAddRun) = this_exitflag;
        add_gpinfo(iAddRun) = this_gpinfo;
    end
    max_x = vertcat(max_x, add_max_x);
    nll = horzcat(nll, add_nll);
    exitflag = horzcat(exitflag, add_exitflag);
    gpinfo = horzcat(gpinfo, add_gpinfo);
end

%% Re-evaluate all runs..

parfor iEval=1:length(nll)
    chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-eval-%d.mat', input_id, iEval));
    if exist(chkpt, 'file')
        fprintf('fitModelBADS :: loading eval %d\n', iEval);
        ld = load(chkpt);
        this_optim_params = ld.this_optim_params;
    else
        fprintf('fitModelBADS :: starting eval %d\n', iEval);
        this_optim_params = Fitting.setParamsFields(base_params, fields, max_x(iEval, :));
        [~, ll, ll_var] = eval_fun(this_optim_params, distribs, signals, choices, final_eval_args{:});
        for iP=1:length(this_optim_params)
            this_optim_params(iP).fit_fields = fields;
            this_optim_params(iP).fit_method = 'BADS';
            this_optim_params(iP).ll = ll;
            this_optim_params(iP).ll_var = ll_var;
        end
        fprintf('fitModelBADS :: saving eval %d\n', iEval);
        saveWrapper(chkpt, struct('this_optim_params', this_optim_params));
    end
    optim_results{iEval} = this_optim_params;
end
end

function saveWrapper(filename, str)
save(filename, '-struct', 'str');
end