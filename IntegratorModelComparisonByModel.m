function [bestfits, model_names] = IntegratorModelComparisonByModel(which_base, which_phase, reference, include_null)

switch which_base
    case 'is'
        % Set base parameters for importance sampling model
        true_params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', .1, 'temperature', 0.1, 'lapse', .005, 'trials', 1600, 'updates', 5, 'samples', 5);
        % Set task parameters... the following are based on earlier runs of Model.getThresholdPoints
        % to guarantee ~70% performance, but for consistency across machines the parameters are
        % hard-coded here (whatever parameters are used here enter into the 'uid' for caching the
        % model-fitting results, which may be precomputed remotely). If the 'true_params' above are
        % changed, then these values will need to be changed as well.
        switch which_phase
            case 'lshc'
                true_params.category_info = 0.906;
                true_params.p_match = 0.906;
                true_params.sensory_info = 0.650;
                true_params.var_s = 13.42;
            case 'hslc'
                true_params.category_info = 0.628;
                true_params.p_match = 0.628;
                true_params.sensory_info = .900;
                true_params.var_s = 1.22;
        end
    case 'vb'
        % Set base parameters for importance sampling model
        true_params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'temperature', 0.1, 'lapse', .005, 'trials', 1600, 'updates', 5, 'step_size', 0.05);
        % See comment in 'is' block above
        switch which_phase
            case 'lshc'
                true_params.category_info = 0.908;
                true_params.p_match = 0.908;
                true_params.sensory_info = 0.655;
                true_params.var_s = 12.61;
            case 'hslc'
                true_params.category_info = 0.691;
                true_params.p_match = 0.691;
                true_params.sensory_info = .917;
                true_params.var_s = 1.04;
        end
    case 'itb'
        % Set base parameters for integration-to-bound model
        true_params = Model.newModelParams('model', 'itb', 'var_x', 0.1, 'gamma', .1, 'temperature', 0.4, 'lapse', .005, 'trials', 1600, 'updates', 1, 'bound', 1);
        % See comment in 'is' block above
        switch which_phase
            case 'lshc'
                true_params.category_info = 0.913;
                true_params.p_match = 0.913;
                true_params.sensory_info = 0.674;
                true_params.var_s = 9.89;
            case 'hslc'
                true_params.category_info = 0.650;
                true_params.p_match = 0.650;
                true_params.sensory_info = 0.906;
                true_params.var_s = 1.15;
        end
end

% Make results repeatable
true_params.seed = 24781390;

%% Plot PK of the true model and print some diagnostics

Model.plotPK(true_params, [1 0 0]);
title('PK of true model');

[data, true_cat] = Model.genDataWithParams(true_params);
res = Model.runVectorized(true_params, data);
fprintf('=== True Model [%s] Stats ===\n', upper(true_params.model));
fprintf('\t%d trials x %d frames\n', true_params.trials, true_params.frames);
fprintf('\tsensory info = %.2f\n', true_params.sensory_info);
fprintf('\tcategory info = %.2f\n', true_params.category_info);
fprintf('\tpercent correct = %.1f%%\n', 100*mean(res.choices == true_cat));

%% Fit various 'integrator' models

seed = 2487429;
uid = Model.getModelStringID(true_params);
sigs = Model.logLikelihoodOdds(true_params, data);
[bestfits, ll_train, ll_test, model_names, ~] = IntegratorModelComparison(sigs, res.choices, 50, uid, seed);

if nargin >= 4 && include_null
    % For reference, compute the log likelihood under the null model
    ll_null = -log(2);
    ll_train(end+1, :) = ll_null;
    ll_test(end+1, :) = ll_test;
    model_names{end+1} = 'null';
end

ref_model = strcmpi(model_names, reference);
nModels = length(model_names);

%% Compute LL _differences_ across reruns

% Each model saw the same data every time so their LLs are correlated across runs. The per-datum
% difference is more informative than the difference of overall performance. First column is LL
% difference on training set:
ll_diff = ll_train - ll_train(ref_model,:);
est_ll_diff(:,1) = mean(ll_diff, 2);
sem_ll_diff(:,1) = std(ll_diff, [], 2) ./ sqrt(size(ll_diff, 2));

% Second column is LL difference on test set:
ll_diff = ll_test - ll_test(ref_model,:);
est_ll_diff(:,2) = mean(ll_diff, 2);
sem_ll_diff(:,3) = std(ll_diff, [], 2) ./ sqrt(size(ll_diff, 2));

%% Plot result
hold on;
bar(1:nModels, est_ll_diff);
errorbar((1:nModels)-.15, est_ll_diff(:,1), sem_ll_diff(:,1), 'ok');
errorbar((1:nModels)+.15, est_ll_diff(:,2), sem_ll_diff(:,2), 'ok');
legend({'train', 'test'}, 'Location', 'Northwest');
set(gca, 'XTick', 1:nModels, 'XTickLabel', model_names);
grid on;
ylabel(['\Delta LL from ' reference ' fit']);
title(sprintf('Model fits to %s [%s]', upper(which_base), upper(which_phase)));
end