function [bestfit, ll_train, ll_test, model_names, model_fields] = IntegratorModelComparison(signals, choices, nFolds, prefix, seed)

model_fields = {
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'bound'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'gamma'}
    {'prior_C', 'log_temperature', 'log_lapse_1', 'log_lapse_2', 'bound', 'gamma'}
};

model_names = {'ideal', 'itb', 'gamma', 'itb-gamma'};

% Set RNG for repeatable folds
if ~exist('seed', 'var') || isempty(seed), seed = 793412; end
rng(seed, 'twister');
nTrials = length(choices);
idxShuffle = randperm(nTrials);

idxFoldStart = round(linspace(1, nTrials+1, nFolds+1));
idxFoldEnd = idxFoldStart(2:end)-1;

use_cache =  nargin >=4 && ~isempty(prefix);
for iModel=length(model_fields):-1:1
    parfor iFold=1:nFolds
        test_idx = idxShuffle(idxFoldStart(iFold):idxFoldEnd(iFold));
        train_idx = setdiff(1:nTrials, test_idx);
        if use_cache
            [bestfit(iModel, iFold), ~, ll_train(iModel, iFold), ll_test(iModel, iFold)] = LoadOrRun(@Fitting.fitIntegratorModel, ...
                {signals(train_idx, :), choices(train_idx), model_fields{iModel}, [], [], signals(test_idx, :), choices(test_idx)}, ...
                fullfile('../Precomputed', ['integrator-xval-' prefix '-' model_names{iModel} '-fold=' num2str(iFold) 'of' num2str(nFolds) '-seed=' num2str(seed) '.mat']));
        else
            [bestfit(iModel, iFold), ~, ll_train(iModel, iFold), ll_test(iModel, iFold)] = ...
                Fitting.fitIntegratorModel(signals(train_idx, :), choices(train_idx), model_fields{iModel}, [], [], signals(test_idx, :), choices(test_idx));
        end
    end
end
end