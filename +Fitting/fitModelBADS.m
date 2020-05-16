function [optim_results, exitflag, gpinfo, metadata] = fitModelBADS(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELBADS Search for maximum likelhood model using Bayesian Adaptive Direct Searach
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
chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-grid.mat', input_id));

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
% See @Fitting.fitModelBADS

% To gracefully save/load from partially computed results, we initialize to NaN values, (maybe) load
% from a file, then only evaluate those rows that are still NaN, with periodic checkpointing.
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
needs_eval = find(isnan(loglike));
for ii=1:length(needs_eval)
    iVal = needs_eval(ii);
    [logpdf(iVal), loglike(iVal), variance(iVal)] = ...
        eval_fun(Fitting.setParamsFields(base_params, fields, grid_points(iVal, :)), distribs, signals, choices, eval_args{:});
    
    % Periodically print diagnostics and save things
    if mod(ii, 100) == 0 || ii == length(needs_eval)
        telapse = toc(tstart);
        fprintf('%s :: eval %d of %d\tETA=%.1fs\n', chkpt, ii, length(needs_eval), (length(needs_eval)-ii)*telapse/ii);
        save(chkpt, 'grid_points', 'logpdf', 'loglike', 'variance');
    end
end

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
end

nRuns = 15;
for iRun=nRuns:-1:1
    % Use the top 'nRuns' points from the grid for initialization
    chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-search-%d.mat', input_id, iRun));
    if exist(chkpt, 'file')
        ld = load(chkpt);
        max_x = ld.max_x;
        nll = ld.nll;
        exitflag = ld.exitflag;
        gpinfo = ld.gpinfo;
        fprintf('fitModelBADS :: loading search iteration %d\n', iRun);
    else
        fprintf('fitModelBADS :: starting search iteration %d\n', iRun);
        [max_x(iRun,:), nll(iRun), exitflag(iRun), ~, ~, gpinfo(iRun)] = bads(...
            @(x) -eval_fun(Fitting.setParamsFields(base_params, fields, x), distribs, signals, choices, eval_args{:}), ...
            grid_points(sort_ll(iRun), :), LB, UB, PLB,PUB, [], opts);
        save(chkpt, 'max_x', 'nll', 'exitflag', 'gpinfo');
    end
end

%% Store all of the (valid) runs..

valid = exitflag >= 0;
max_x = max_x(valid, :);
nll = nll(valid);

for iRun=length(nll):-1:1
    chkpt = fullfile('sample-checkpoints', sprintf('%x-bads-eval-%d.mat', input_id, iRun));
    if exist(chkpt, 'file')
        ld = load(chkpt);
        optim_results = ld.optim_results;
        fprintf('fitModelBADS :: loading eval %d\n', iRun);
    else
        optim_results{iRun} = Fitting.setParamsFields(base_params, fields, max_x(iRun, :));
        [~, ll, ll_var] = eval_fun(optim_results{iRun}, distribs, signals, choices, final_eval_args{:});
        for iP=1:length(optim_results{iRun})
            optim_results{iRun}(iP).fit_fields = fields;
            optim_results{iRun}(iP).fit_method = 'BADS';
            optim_results{iRun}(iP).ll = ll;
            optim_results{iRun}(iP).ll_var = ll_var;
        end
        fprintf('fitModelBADS :: saving eval %d\n', iRun);
        save(chkpt, 'optim_results');
    end
end
end