function [optim_results, samples, sample_scores, metadata] = fitModelMH(base_params, signals, choices, distribs, metadata)
%FITTING.FITMODELMH Strategy for fitting inference models by MH sampling followed by selecting the
%best sample as an approximation to the MAP/MLE. Stochastic models are handled as a special case.

fields = fieldnames(distribs);

% Get a UID for the model being fit, signals, and choices so that we can restart from checkpoints
if iscell(signals)
    allsigs = vertcat(signals{:});
    allchoices = vertcat(choices{:});
    input_id = string2hash([base_params(1).model, strjoin(fieldnames(distribs)), num2str([allsigs(:)' allchoices'])]);
    nTrials = length(allchoices);
else
    input_id = string2hash([base_params(1).model, strjoin(fieldnames(distribs)), num2str([signals(:)' choices'])]);
    nTrials = length(choices);
end
chkpt = fullfile('sample-checkpoints', sprintf('%x', input_id));

stoch = Model.isStochastic(base_params);
if stoch
    % Fitting strategy 1 (stochastic model): (A) MH sample points, resulting in over-dispersed
    % samples due to internal stochasticity (B) roughly score each sample with IBS, then (C) more
    % carefully search for the best one.
    
    % Step 1A: sample. Note that stochastic models with finite repetitions will result in posteriors
    % that are wider than the true posterior because each LH evaluation is corrupted by noise
    disp('fitModelMH: stochastic step A');
    [samples, ~, sample_scores] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 5, 50, 10, chkpt);
    
    % Step 1B: Get a rough estimate of the log probabilities at each (unique) sample
    disp('fitModelMH: stochastic step B');
    stretch = 3;
    [uSamples, ~, idxExpand] = unique(samples, 'rows');
    for iSamp=size(uSamples, 1):-1:1
        est_ll(iSamp) = mean(sample_scores.loglike(idxExpand == iSamp));
        est_lp(iSamp) = mean(sample_scores.logpdf(idxExpand == iSamp));
        est_var(iSamp) = mean(sample_scores.variance(idxExpand == iSamp)) / sum(idxExpand == iSamp);
    end
    upper_conf = est_ll + stretch*sqrt(est_var);
    lower_conf = est_ll - stretch*sqrt(est_var);
    
    % Step 1C: Search for the best sample (maximum likelihood)
    disp('fitModelMH: stochastic step C');
    
    refine_file = [chkpt '-refine.mat'];
    if exist(refine_file, 'file')
        ld = load(refine_file);
        samples = ld.samples;
        sample_scores = ld.sample_scores;
        est_ll = sample_scores.loglike;
        est_lp = sample_scores.logpdf;
        est_var = sample_scores.variance;
    else
        allidx = 1:size(uSamples,1);
        [~, ibest] = max(est_ll);
        
        max_eval_per_point = round(sqrt(nTrials));
        n_evals_per = ones(size(upper_conf));
        
        % Loop until our *lower* bound on the best point is better than the *upper* bound on all other
        % points
        while lower_conf(ibest) < max(upper_conf(setdiff(allidx, ibest)))
            % Determine which point to evaluate further: choose the point with the highest upper
            % bound as long as it is within range of the current best point's lower bound, excluding
            % those that have already been evaluated 'max_eval_per_point' times.
            eligible = upper_conf > lower_conf(ibest) & n_evals_per < max_eval_per_point;
            if ~any(eligible), break; end
            tmp_uc = upper_conf;
            tmp_uc(~eligible) = -inf;
            [~, ieval] = max(tmp_uc);
            
            % DEBUG FIGURE
            cla; hold on;
            errorbar(1:size(uSamples,1), est_ll, 2*sqrt(est_var));
            [~,ibest] = max(est_ll);
            errorbar(ibest, est_ll(ibest), 2*sqrt(est_var(ibest)), '-r');
            plot(ieval, est_ll(ieval), 'og');
            drawnow;

            % Compute a better estimate of the log prob for this point by using a larger number of
            % repeats.
            [new_est_lp, new_est_ll, new_est_var] = Fitting.choiceModelLogProbIBS(...
                Fitting.setParamsFields(base_params, fields, uSamples(ieval, :)), distribs, signals, choices, [], 5);
            
            % Update stuff
            n_evals_per(ieval) = n_evals_per(ieval) + 1;
            
            est_ll(ieval) = (est_ll(ieval)*new_est_var + new_est_ll*(est_var(ieval))) / (est_var(ieval) + new_est_var);
            est_lp(ieval) = (est_lp(ieval)*new_est_var + new_est_lp*(est_var(ieval))) / (est_var(ieval) + new_est_var);
            est_var(ieval) = 1./(1./est_var(ieval) + 1./new_est_var);
            
            upper_conf(ieval) = est_ll(ieval) + stretch*sqrt(est_var(ieval));
            lower_conf(ieval) = est_ll(ieval) - stretch*sqrt(est_var(ieval));
        end
        
        % Store refined estimates back into 'sample_scores'
        sample_scores.loglike = est_ll(idxExpand);
        sample_scores.logpdf = est_lp(idxExpand);
        sample_scores.variance = est_var(idxExpand);
        
        save(refine_file, 'samples', 'sample_scores');
    end
    
    gp_args = {'Sigma', sqrt(mean(est_var)), 'SigmaLowerBound', min(sqrt(est_var)), 'ConstantSigma', false};
else
    % Fitting strategy 2 (deterministic model): (A) sample the posterior then (B) pick the best
    % sample.
    
    % Step 2A: sample
    disp('fitModelMH: deterministic step A');
    [samples, ~, sample_scores] = Fitting.sampleModelMH(signals, choices, base_params, 2000, distribs, 0, 50, 10, chkpt);
    
    % Step 2B: prep GP below by combining all evals for redundant samples
    [uSamples, ~, idxExpand] = unique(samples, 'rows');
    for iSamp=size(uSamples, 1):-1:1
        est_ll(iSamp) = mean(sample_scores.loglike(idxExpand == iSamp));
        est_lp(iSamp) = mean(sample_scores.logpdf(idxExpand == iSamp));
        est_var(iSamp) = mean(sample_scores.variance(idxExpand == iSamp)) / sum(idxExpand == iSamp);
    end
    
    gp_args = {'Sigma', .01, 'ConstantSigma', true};
end
        
% Penultimate step: find maximum likelihood and MAP sample
[~, best_ll_idx] = max(sample_scores.loglike);
optim_results.mle_params = Fitting.setParamsFields(base_params, fields, samples(best_ll_idx, :));

[~, best_lp_idx] = max(sample_scores.logpdf);
optim_results.map_params = Fitting.setParamsFields(base_params, fields, samples(best_lp_idx, :));

%% Finally, search LL landscape with GP, starting from the best point

% Set initial length scales as 1/20th of the range of plausible bounds for each parameter
init_scale = cellfun(@(f) (distribs.(f).pub-distribs.(f).plb)/20, fields)';

% Fit Gaussian process to all log likelihood evaluations
gp_ll = fitrgp(uSamples, est_ll, 'KernelFunction', 'ardsquaredexponential', ...
    'KernelParameters', [init_scale std(est_ll)], gp_args{:});
% Search the GP fit for a better maximum
gp_mle = fmincon(@(x) -gp_ll.predict(x), samples(best_ll_idx, :), [], [], [], [], ...
    cellfun(@(f) distribs.(f).lb, fields), cellfun(@(f) distribs.(f).ub, fields));
optim_results.gp_mle_params = Fitting.setParamsFields(base_params, fields, gp_mle);

% Repeat the above for the MAP
gp_lp = fitrgp(uSamples, est_lp, 'KernelFunction', 'ardsquaredexponential', ...
    'KernelParameters', [init_scale std(est_lp)], gp_args{:});
% Search the GP fit for a better maximum
gp_map = fmincon(@(x) -gp_lp.predict(x), samples(best_lp_idx, :), [], [], [], [], ...
    cellfun(@(f) distribs.(f).lb, fields), cellfun(@(f) distribs.(f).ub, fields));
optim_results.gp_map_params = Fitting.setParamsFields(base_params, fields, gp_map);

end