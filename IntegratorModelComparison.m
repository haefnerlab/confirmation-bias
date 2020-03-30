function [bestfit, ll_train, ll_test, n_train, model_names, model_fields] = IntegratorModelComparison(signals, choices, fracTrain, nRepeats, prefix)

model_fields = {
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'bound'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'gamma'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'bound', 'gamma'}
};

model_names = {'ideal', 'itb', 'gamma', 'itb_gamma'};

use_cache =  nargin >= 5 && ~isempty(prefix);
parfor iRep=1:nRepeats
    if use_cache
        [bestfit(:, iRep), ll_train(:, iRep), ll_test(:, iRep), n_train(iRep)] = LoadOrRun(@fit_wrapper, ...
            {signals, choices, model_fields, fracTrain}, ...
            fullfile('../Precomputed', ['integrator-xval-' prefix '-[' strjoin(model_names, '-') ']-run=' num2str(iRep) '.mat']));
    else
        [bestfit(:, iRep), ll_train(:, iRep), ll_test(:, iRep), n_train(iRep)] = fit_wrapper(signals, choices, model_fields, fracTrain);
    end
end
end

function [bestfits, ll_trains, ll_tests, nTrain] = fit_wrapper(signals, choices, model_fields, fracTrain)
% Run all models on the same random reshuffling of the data. If 'fracTrain' is empty, we randomly
% decide how many to train on (approximating leave-p-out averaged over all p)
if isempty(fracTrain), fracTrain = rand(); end
nTrials = length(choices);
nTrain = min(round(fracTrain * nTrials), nTrials-1);
idx_shuffle = randperm(nTrials);
train_idx = idx_shuffle(1:nTrain);
test_idx = idx_shuffle(nTrain+1:end);

for iModel=length(model_fields):-1:1
    [bestfits(iModel, 1), ~, ll_trains(iModel, 1), ll_tests(iModel, 1)] = ...
        Fitting.fitIntegratorModel(signals(train_idx, :), choices(train_idx), model_fields{iModel}, [], [], signals(test_idx, :), choices(test_idx));
end
end