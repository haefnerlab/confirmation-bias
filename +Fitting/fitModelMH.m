function [map_params, samples, scores, metadata] = fitModelMHMAP(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELMHMAP Strategy for fitting inference models by MH sampling followed by selecting
%the best sample as an approximation to the MAP. Stochastic models are handled as a special case.

fields = fieldnames(distribs);

stoch = Model.isStochastic(base_params);
if stoch
    % Fitting strategy 1 (stochastic model): (A) MH sample points, (B) roughly score them with IBS,
    % then (C) more carefully search for the best one.
    
    % Step 1A: sample. Note that stochastic models with finite 'inner loop' simulations will result
    % in posteriors that are wider than the true posterior because each LH evaluation is corrupted
    % by noise
    disp('fitModelMHMAP: stochastic step A');
    input_id = string2hash([Model.getModelStringID(base_params) num2str([signals(:)' choices'])]);
    chkpt = fullfile('sample-checkpoints', sprintf('%x', input_id));
    samples = Fitting.fitChoicesMH(signals, choices, base_params, 20000, distribs, 10, 0, 1, chkpt);
    
    % Step 1B: Get a rough estimate of the posterior probability at each (unique) sample
    disp('fitModelMHMAP: stochastic step B');
    [uSamples, ~, idxExpand] = unique(samples, 'rows');
    for iSamp=size(uSamples,1):-1:1
        [est_log_post(iSamp), ~, est_var(iSamp)] = Fitting.choiceModelLogProbIBS(Fitting.setParamsFields(base_params, fields, uSamples(iSamp, :)), ...
            distribs, signals, choices);
    end
    upper_conf = est_log_post+3*sqrt(est_var);
    
    % Step 1C: Search for the best sample
    disp('fitModelMHMAP: stochastic step C');
    [srt_upper_conf, isrt] = sort(upper_conf, 'descend');
    uSamples = uSamples(isrt, :);
    
    precise_logprob = -inf(size(est_log_post));
    best_idx = 1;
    for iSamp=1:size(uSamples,1)
        % Break when upper confidence region of all remaining samples is below the best estimate
        if precise_logprob(best_idx) > srt_upper_conf(iSamp)
            break;
        end
        
        % Compute a better estimate of the log prob for this point by using a larger number of 'inner
        % loop' simulations per point.
        precise_logprob(iSamp) = Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, uSamples(iSamp, :)), ...
            distribs, signals, choices, 150);
        
        if precise_logprob(iSamp) > precise_logprob(best_idx)
            best_idx = iSamp;
        end
    end
    
    % DEBUG FIGURE: plot log prob bounds and actual values in order they were checked
    figure; hold on;
    errorbar(1:size(uSamples,1), est_log_post(isrt), sqrt(est_var(isrt)));
    errorbar(1:size(uSamples,1), est_log_post(isrt), 4*sqrt(est_var(isrt)));
    plot(precise_logprob);
    
    map_params = Fitting.setParamsFields(base_params, fields, uSamples(best_idx, :));
    invsrt = zeros(size(isrt));
    invsrt(isrt) = 1:length(isrt);
    precise_logprob = precise_logprob(invsrt);
    scores = precise_logprob(idxExpand);
else
    % Fitting strategy 2 (deterministic model): (A) sample the posterior then (B) pick the best
    % sample.
    
    % Step 2A: sample
    disp('fitModelMHMAP: deterministic step A');
    samples = Fitting.fitChoicesMH(signals, choices, base_params, 2000, distribs, 1);
    
    % Step 2B: compute log posterior score of each sample
    disp('fitModelMHMAP: deterministic step B');
    [uSamples, ~, idxExpand] = unique(samples, 'rows');
    for iSamp=size(uSamples,1):-1:1
        log_post(iSamp) = Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, uSamples(iSamp, :)), ...
            distribs, signals, choices, 1);
    end
    
    [~, best_idx] = max(log_post);
    map_params = Fitting.setParamsFields(base_params, fields, uSamples(best_idx, :));
    scores = log_post(idxExpand);
end


end