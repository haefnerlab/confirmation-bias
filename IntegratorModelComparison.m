function [ll_train, ll_test, n_train, model_names, model_fields] = IntegratorModelComparison(signals, choices, fracTrainRange, nSplits, prefix)

model_fields = {
    {'prior_C', 'log_temperature', 'log_lapse1', 'log_lapse2'}
    {'prior_C', 'log_temperature', 'log_lapse1', 'log_lapse2', 'bound'}
    {'prior_C', 'log_temperature', 'log_lapse1', 'log_lapse2', 'gamma'}
    {'prior_C', 'log_temperature', 'log_lapse1', 'log_lapse2', 'bound', 'gamma'}
};

model_names = {'ideal', 'itb', 'gamma', 'itb_gamma'};

use_cache =  nargin >= 5 && ~isempty(prefix);
parfor iSplit=1:nSplits
    if use_cache
        [~, ll_train(:, iSplit), ll_test(:, iSplit), n_train(iSplit)] = ...
            LoadOrRun(@run_split_wrapper, {signals, choices, model_fields, fracTrainRange}, ...
            fullfile('../Precomputed', ['integrator-xval-' prefix '-[' strjoin(model_names, '-') ']-run=' num2str(iSplit) '.mat']));
    else
        [~, ll_train(:, iSplit), ll_test(:, iSplit), n_train(iSplit)] = ...
            run_split_wrapper(signals, choices, model_fields, fracTrainRange);
    end
end
end

function [samples, ll_trains, ll_tests, nTrain] = run_split_wrapper(signals, choices, model_fields, fracTrainRange)
% First, decide how many data points will be in the train vs test set for this run. We allow
% 'fracTrain' to be either a scalar to set p to a single value (as in classic leave-p-out for a
% fixed 'p') or a range defining [plow phi], e.g. [0.1 1.0] to set minimum and maximum % of data
% left out. The 'true' model score is with fracTrain=[0 1] to get all data points, but setting a
% minimum like [0.1 1.0] helps stabilize the result [1].
%
% [1] Fong, E., & Holmes, C. (2019). On the marginal likelihood and cross-validation. Retrieved from
%     http://arxiv.org/abs/1905.08737
thisFracTrain = min(fracTrainRange) + rand()*(max(fracTrainRange)-min(fracTrainRange));

nTrials = length(choices);
nTrain = min(round(thisFracTrain * nTrials), nTrials-1);
idx_shuffle = randperm(nTrials);
train_idx = idx_shuffle(1:nTrain);
test_idx = idx_shuffle(nTrain+1:end);

for iModel=length(model_fields):-1:1
    [samples{iModel}, ~, base_params] = Fitting.sampleIntegratorModel(signals(train_idx,:), ...
        choices(train_idx), 1000, model_fields{iModel}, true, [], 'burnin', 500, 'thin', 25);
    
    % Evaluate [samples x data points] log likelihood values
    for s=size(samples{iModel}):-1:1
        run_results = Model.runIntegratorModel(...
            Fitting.setParamsFields(base_params, model_fields{iModel}, samples{iModel}(s,:)), signals);
        likelihoods(s, :) = run_results.prob_choice(choices == +1) + run_results.prob_choice(choices ~= +1);
    end
    
    % Average likelihoods across samples per data point to estimate log likelihood. Note: slight
    % bias in this estimate
    log_likelihoods = log(mean(likelihoods, 1));
    
    ll_trains(iModel) = mean(log_likelihoods(train_idx));
    ll_tests(iModel) = mean(log_likelihoods(test_idx));
end
end