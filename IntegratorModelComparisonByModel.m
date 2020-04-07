function [true_params, model_names] = IntegratorModelComparisonByModel(which_base, which_phase, group_by, reference, add_null)

if nargin < 3, group_by = 'traintest'; end % or 'model'
if nargin < 4, reference = 'ideal'; end
if nargin < 5, add_null = false; end

base_params = Model.newModelParams('var_x', 0.1, 'gamma', .1, 'temperature', 0.1, 'lapse', .005, 'trials', 1600);
base_params.save_dir = 'tmp';
switch which_base
    case 'is'
        % Set parameters for importance sampling model
        true_params = base_params;
        true_params.model = 'is';
        true_params.updates = 5;
        true_params.samples = 5;
    case 'vb'
        % Set parameters for VB model
        true_params = base_params;
        true_params.model = 'vb-czx';
        true_params.updates = 5;
        true_params.step_size = 0.05;
    case 'itb'
        % Set base parameters for integration-to-bound model
        true_params = base_params;
        true_params.model = 'itb';
        true_params.updates = 1;
        true_params.bound = 1.2;
        true_params.noise = 0.35;
        true_params.gammafun = @(ci,si) (1-ci);
end

sens_cat_pts = Model.getThresholdPoints(0.51:0.02:0.99, true_params, 0.7, 5);

switch which_phase
    case 'lshc'
        si = sens_cat_pts(1,1);
        ci = sens_cat_pts(1,2);
    case 'hslc'
        si = sens_cat_pts(end,1);
        ci = sens_cat_pts(end,2);
end

true_params = Model.setCategorySensoryInfo(true_params, ci, si);

% Make results repeatable
true_params.seed = 977116605;

%% Print some diagnostics about the true model

[data, true_cat] = Model.genDataWithParams(true_params);
res = Model.runVectorized(true_params, data);
fprintf('=== True Model [%s] Stats ===\n', upper(true_params.model));
fprintf('\t%d trials x %d frames\n', true_params.trials, true_params.frames);
fprintf('\tsensory info = %.2f\n', true_params.sensory_info);
fprintf('\tcategory info = %.2f\n', true_params.category_info);
fprintf('\tpercent correct = %.1f%%\n', 100*mean(res.choices == true_cat));
[~,~,abb] = Model.plotPK(true_params, 'exp', {}, [], false);
fprintf('\tPK beta fit = %.1f\n', abb(2));

%% Fit various 'integrator' models

uid = Model.getModelStringID(true_params);
sigs = Model.logLikelihoodOdds(true_params, data);
[ll_train, ll_test, ~, model_names] = IntegratorModelComparison(sigs, res.choices, [0.1 1], 500, uid);

if add_null
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
sem_ll_diff(:,2) = std(ll_diff, [], 2) ./ sqrt(size(ll_diff, 2));

%% Plot result
hold on;
if startsWith(group_by, 'm')
    h = bar(est_ll_diff); drawnow;
    errorbar(h(1).XData+h(1).XOffset, est_ll_diff(:,1), sem_ll_diff(:,1), 'ok');
    errorbar(h(2).XData+h(2).XOffset, est_ll_diff(:,2), sem_ll_diff(:,2), 'ok');
    legend({'train', 'test'}, 'Location', 'Northwest');
    set(gca, 'XTick', 1:nModels, 'XTickLabel', model_names);
elseif startsWith(group_by, 't')
    h = bar(est_ll_diff'); drawnow;
    for m=1:length(model_names)
        errorbar(h(m).XData+h(m).XOffset, est_ll_diff(m,:), sem_ll_diff(m,:), 'ok');
    end
    legend(model_names, 'Location', 'Northwest');
    set(gca, 'XTick', 1:2, 'XTickLabel', {'train', 'test'});
end
grid on;
ylabel(['\Delta LL from ' reference ' fit']);
title(sprintf('Model fits to %s [%s]', upper(which_base), upper(which_phase)));
end