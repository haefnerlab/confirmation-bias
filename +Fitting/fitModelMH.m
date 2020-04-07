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
    if iscell(signals)
        allsigs = vertcat(signals{:});
        allchoices = vertcat(choices{:});
        input_id = string2hash([Model.getModelStringID(base_params(1)) num2str([allsigs(:)' allchoices'])]);
    else
        input_id = string2hash([Model.getModelStringID(base_params) num2str([signals(:)' choices'])]);
    end
    chkpt = fullfile('sample-checkpoints', sprintf('%x', input_id));
    [samples, ~, sample_info] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 5, 50, 10, chkpt);
    
    % Step 1B: Get a rough estimate of the log probabilities at each (unique) sample
    disp('fitModelMH: stochastic step B');
    [uSamples, ~, idxExpand] = unique(samples, 'rows');
    for iSamp=size(uSamples, 1):-1:1
        est_ll(iSamp) = mean(sample_info.loglike(idxExpand == iSamp));
        est_lp(iSamp) = mean(sample_info.logpdf(idxExpand == iSamp));
        est_var(iSamp) = mean(sample_info.variance(idxExpand == iSamp)) / sum(idxExpand == iSamp);
    end
    upper_conf = est_ll + 2*sqrt(est_var);
    lower_conf = est_ll - 2*sqrt(est_var);
    
    % Step 1C: Search for the best sample (maximum likelihood)
    disp('fitModelMH: stochastic step C');
    
    refine_file = [chkpt '-refine.mat'];
    if exist(refine_file, 'file')
        ld = load(refine_file);
        samples = ld.samples;
        sample_scores = ld.sample_scores;
    else
        allidx = 1:size(uSamples,1);
        [~, imax] = max(upper_conf);
        
        maxeval = 10000;
        ieval = 1;
        % Loop until our *lower* bound on the best point is better than the *upper* bound on all other
        % points
        while lower_conf(imax) < max(upper_conf(setdiff(allidx, imax))) && ieval < maxeval
            % % DEBUG FIGURE
            % cla; hold on;
            % errorbar(1:size(uSamples,1), est_ll, 2*sqrt(est_var));
            % [~,ibest] = max(est_ll);
            % errorbar(ibest, est_ll(ibest), 2*sqrt(est_var(ibest)), '-r');
            % plot(imax, est_ll(imax), 'og');
            % drawnow;

            % Compute a better estimate of the log prob for this point by using a larger number of
            % repeats.
            [new_est_lp, new_est_ll, new_est_var] = Fitting.choiceModelLogProbIBS(...
                Fitting.setParamsFields(base_params, fields, uSamples(imax, :)), distribs, signals, choices, [], 50);
            
            est_ll(imax) = (est_ll(imax)*new_est_var + new_est_ll*(est_var(imax))) / (est_var(imax) + new_est_var);
            est_lp(imax) = (est_lp(imax)*new_est_var + new_est_lp*(est_var(imax))) / (est_var(imax) + new_est_var);
            est_var(imax) = 1./(1./est_var(imax) + 1./new_est_var);
            
            upper_conf(imax) = est_ll(imax) + 3*sqrt(est_var(imax));
            lower_conf(imax) = est_ll(imax) - 3*sqrt(est_var(imax));
            
            [~, imax] = max(upper_conf);
            
            ieval = ieval + 1;
        end
        
        % Store refined estimates back into 'sample_scores'
        sample_scores.loglike = est_ll(idxExpand);
        sample_scores.logpdf = est_lp(idxExpand);
        sample_scores.variance = est_var(idxExpand);
        
        save(refine_file, 'samples', 'sample_scores');
    end
else
    % Fitting strategy 2 (deterministic model): (A) sample the posterior then (B) pick the best
    % sample.
    
    % Step 2A: sample
    disp('fitModelMH: deterministic step A');
    [samples, ~, sample_scores] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 0, 50, 10);
end
        
% Penultimate step: find maximum likelihood and MAP sample
[~, best_ll_idx] = max(sample_scores.loglike);
mle_params = Fitting.setParamsFields(base_params, fields, samples(best_ll_idx, :));

[~, best_lp_idx] = max(sample_scores.logpdf);
map_params = Fitting.setParamsFields(base_params, fields, samples(best_lp_idx, :));

%% Finally, search LL landscape with GP, starting from the best point


end