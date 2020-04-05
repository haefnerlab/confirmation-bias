function [samples, mle_params, map_params, sample_scores, metadata] = fitModelMH(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELMH Strategy for fitting inference models by MH sampling followed by selecting the
%best sample as an approximation to the MAP/MLE. Stochastic models are handled as a special case.

fields = fieldnames(distribs);

stoch = Model.isStochastic(base_params);
if stoch
    % Fitting strategy 1 (stochastic model): (A) MH sample points, resulting in over-dispersed
    % samples due to internal stochasticity (B) roughly score each sample with IBS, then (C) more
    % carefully search for the best one.
    
    % Step 1A: sample. Note that stochastic models with finite repetitions will result in posteriors
    % that are wider than the true posterior because each LH evaluation is corrupted by noise
    disp('fitModelMH: stochastic step A');
    input_id = string2hash([Model.getModelStringID(base_params) num2str([signals(:)' choices'])]);
    chkpt = fullfile('sample-checkpoints', sprintf('%x', input_id));
    [samples, ~, sample_info] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 5, 50, 10, chkpt);
    
    % Step 1B: Get a rough estimate of the posterior probability at each (unique) sample
    disp('fitModelMH: stochastic step B');
    [uSamples, ~, idxExpand] = unique(samples, 'rows');
    for iSamp=size(uSamples, 1):-1:1
        est_ll(iSamp) = mean(sample_info.loglike(idxExpand == iSamp));
        est_lp(iSamp) = mean(sample_info.logpdf(idxExpand == iSamp));
        est_var(iSamp) = mean(sample_info.variance(idxExpand == iSamp)) / sum(idxExpand == iSamp);
    end
    upper_conf = est_ll + 3*sqrt(est_var);
    
    % Step 1C: Search for the best sample (maximum likelihood)
    disp('fitModelMH: stochastic step C');
    [srt_upper_conf, isrt] = sort(upper_conf, 'descend');
    uSamples = uSamples(isrt, :);
    est_ll = est_ll(isrt);
    est_lp = est_lp(isrt);
    est_var = est_var(isrt);
    
    best_ll_idx = 1;
    best_lp_idx = 1;
    for iSamp=1:size(uSamples,1)
        % Break when upper confidence region of all remaining samples is below the lower confidence
        % region of the best estimate (LL)
        if est_ll(best_ll_idx)-3*sqrt(est_var(best_idx)) > srt_upper_conf(iSamp)
            break;
        end
        
        % Compute a better estimate of the log prob for this point by using a larger number of
        % repeats.
        [new_est_lp, new_est_ll, new_est_var] = Fitting.choiceModelLogProbIBS(...
            Fitting.setParamsFields(base_params, fields, uSamples(iSamp, :)), distribs, signals, choices, [], 50);
        
        est_ll(iSamp) = (est_ll(iSamp)*new_est_var + new_est_ll*(est_var(iSamp))) / (est_var(iSamp) + new_est_var);
        est_lp(iSamp) = (est_lp(iSamp)*new_est_var + new_est_lp*(est_var(iSamp))) / (est_var(iSamp) + new_est_var);
        est_var(iSamp) = 1./(1./est_var(iamp) + 1./new_est_var);
        
        if est_ll(iSamp) > est_ll(best_ll_idx)
            best_ll_idx = iSamp;
        end
        
        if est_lp(iSamp) > est_ll(best_lp_idx)
            best_lp_idx = iSamp;
        end
    end
    
    % DEBUG FIGURE
    figure; hold on;
    errorbar(1:size(uSamples,1), est_ll(isrt), sqrt(est_var(isrt)));
    errorbar(1:size(uSamples,1), est_ll(isrt), 3*sqrt(est_var(isrt)));
    plot(best_ll_idx*[1 1], ylim, '--r');
    
    % Store best post and best like
    mle_params = Fitting.setParamsFields(base_params, fields, uSamples(best_ll_idx, :));
    map_params = Fitting.setParamsFields(base_params, fields, uSamples(best_lp_idx, :));

    % Un-sort the values and return estimates to [nsamples x ...] shape
    invsrt = zeros(size(isrt));
    invsrt(isrt) = 1:length(isrt);
    est_ll = est_ll(invsrt);
    est_lp = est_lp(invsrt);
    est_var = est_var(invsrt);
    sample_scores.loglike = est_ll(idxExpand);
    sample_scores.logpdf = est_lp(idxExpand);
    sample_scores.variance = est_var(idxExpand);
else
    % Fitting strategy 2 (deterministic model): (A) sample the posterior then (B) pick the best
    % sample.
    
    % Step 2A: sample
    disp('fitModelMH: deterministic step A');
    [samples, ~, sample_scores] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 0, 50, 10);
    
    % Step 2B: find maximum likelihood and MAP sample
    [~, best_ll_idx] = max(sample_scores.loglike);
    mle_params = Fitting.setParamsFields(base_params, fields, samples(best_ll_idx, :));
    
    [~, best_lp_idx] = max(sample_scores.logpdf);
    map_params = Fitting.setParamsFields(base_params, fields, samples(best_lp_idx, :));
end


end