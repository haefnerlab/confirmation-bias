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
chkpt = fullfile('fit-checkpoints', sprintf('%x-bads-grid.mat', input_id));
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

%% Run BADS multiple times

LB = cellfun(@(f) distribs.(f).lb, fields)';
UB = cellfun(@(f) distribs.(f).ub, fields)';
PLB = cellfun(@(f) distribs.(f).plb, fields)';
PUB = cellfun(@(f) distribs.(f).pub, fields)';

% Searches are initialized pseudo-randomly from the existing grid + previous evaluations, plus
% this amount of jittering
jitteramt = (PUB-PLB)/20;

opts = bads('defaults');
opts.NonlinearScaling = 'off';
if Model.isStochastic(base_params)
    opts.UncertaintyHandling = 'on';
    opts.NoiseScale = sqrt(mean(variance(loglike > median(loglike))));
else
    opts.UncertaintyHandling = 'off';
    opts.NoiseScale = 0;
end

% Initialize empty return values (overhead of growing arrays is not a concern)
est_ll = [];
est_ll_var = [];
gpinfo = [];
exitflag = [];
optim_results = {};

% Loop until convergence: use existing set of minima (once we have at least 20) to extrapolate how
% many runs we'd need to be confident in the global minimum (see @Fitting.bootstrapMinimaregret). As
% long as we're < this minimum, run and save another 5 calls to BADS.
maxRuns = 200;
estMinRuns = maxRuns;
while length(est_ll) < min(maxRuns, estMinRuns)
    fprintf('fitModelBADS :: completed %d searches. Extending to %d\n', length(est_ll), estMinRuns);
        
    % Get starting points for the next batch of parfor evaluations by a sort of empirical Thompson
    % sampling + jittering. Searching all points in the QRG plus points that have been evaluated so
    % far. Kind of like Bayesian Optimization of the initial point fed to BADS.
    for iAddRun=5:-1:1
        all_loglike = vertcat(loglike(:), est_ll(:));
        all_loglike_var = vertcat(variance(:), est_ll_var(:));
        [~, ibest] = max(all_loglike+sqrt(all_loglike_var).*randn(size(all_loglike)));
        start_pt(iAddRun,:) = grid_points(ibest, :) + randn(size(jitteramt)).*jitteramt;
        start_pt(iAddRun,:) = max(PLB, min(PUB, start_pt(iAddRun,:)));
    end

    parfor iAddRun=1:5
        iRun = iAddRun+length(est_ll);
        
        %% Search step
        chkpt = fullfile('fit-checkpoints', sprintf('%x-bads-search-%d.mat', input_id, iRun));
        if exist(chkpt, 'file')
            ld = load(chkpt);
            this_optim_x = ld.this_max_x;
            this_exitflag = ld.this_exitflag;
            this_gpinfo = ld.this_gpinfo;
            fprintf('fitModelBADS :: loading search %d\n', iRun);
        else
            fprintf('fitModelBADS :: starting search %d\n', iRun);
            [this_optim_x, this_nll, this_exitflag, ~, ~, this_gpinfo] = bads(...
                @(x) -eval_fun(Fitting.setParamsFields(base_params, fields, x), distribs, signals, choices, eval_args{:}), ...
                start_pt(iAddRun,:), LB, UB, PLB, PUB, [], opts);
            fprintf('fitModelBADS :: saving search %d\n', iRun);
            saveWrapper(chkpt, struct('this_max_x', this_optim_x, 'this_nll', this_nll, ...
                'this_exitflag', this_exitflag, 'this_gpinfo', this_gpinfo));
        end
        add_max_x(iAddRun,:) = this_optim_x;
        add_exitflag(iAddRun) = this_exitflag;
        add_gpinfo(iAddRun) = this_gpinfo;
        
        %% Eval step
        chkpt = fullfile('fit-checkpoints', sprintf('%x-bads-eval-%d.mat', input_id, iRun));
        if exist(chkpt, 'file')
            fprintf('fitModelBADS :: loading eval %d\n', iRun);
            ld = load(chkpt);
            this_optim_params = ld.this_optim_params;
        else
            fprintf('fitModelBADS :: starting eval %d\n', iRun);
            this_optim_params = Fitting.setParamsFields(base_params, fields, add_max_x(iAddRun,:));
            [~, ll, ll_var] = eval_fun(this_optim_params, distribs, signals, choices, final_eval_args{:});
            for iP=1:length(this_optim_params)
                this_optim_params(iP).fit_fields = fields;
                this_optim_params(iP).fit_method = 'BADS';
                this_optim_params(iP).ll = ll;
                this_optim_params(iP).ll_var = ll_var;
            end
            fprintf('fitModelBADS :: saving eval %d\n', iRun);
            saveWrapper(chkpt, struct('this_optim_params', this_optim_params));
        end
        this_optim_results{iAddRun} = this_optim_params;
    end
    % Concatenate previous results with new parfor evals
    grid_points = vertcat(grid_points, add_max_x);
    exitflag = horzcat(exitflag, add_exitflag);
    gpinfo = horzcat(gpinfo, add_gpinfo);
    optim_results = horzcat(optim_results, this_optim_results);
    
    % Grab LL and LL_var fields from each fit and use them to estimate min # runs needed for
    % convergence (as long as we've completed at least 20 so far.)
    est_ll = cellfun(@(fit_para) fit_para(1).ll, optim_results);
    est_ll_var = cellfun(@(fit_para) fit_para(1).ll_var, optim_results);
    if length(est_ll) >= 20
        if Model.isStochastic(base_params)
            % How many runs will it take to be 99% sure that our MLE+/-error contains the max?
            estMinRuns = Fitting.bootstrapMinimaRegret(-est_ll, sqrt(est_ll_var), 0);
        else
            % How many runs will it take to be 99% sure that our MLE is within 1 unit of LL from the
            % global max?
            estMinRuns = Fitting.bootstrapMinimaRegret(-est_ll, [], 1);
        end
    end
end
end

function saveWrapper(filename, str)
save(filename, '-struct', 'str');
end