function ModelComparisonFitToSampling(which_base, which_phase)
ps = 0.51:0.02:0.99;

switch which_base
    case 'is'
        % Run with IS as true model
        params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', .1, 'temperature', 0.1, 'trials', 1600, 'updates', 5, 'samples', 5);
    case 'itb'
        % Run with ITB as true model
        params = Model.newModelParams('model', 'itb', 'var_x', 0.1, 'gamma', .1, 'temperature', 0.4, 'bound', 1, 'trials', 1600, 'updates', 1);
end

[~, sens_cat_pts] = Model.getThresholdPoints(ps, params, .7, 5);
params.seed = 1376289;

switch which_phase
    case 'lshc'
        % Use LSHC as reference, fit other models to it
        true_params = params;
        true_params.sensory_info = sens_cat_pts(1,1);
        true_params.var_s = Model.getEvidenceVariance(sens_cat_pts(1,1));
        true_params.category_info = sens_cat_pts(1,2);
        true_params.p_match = sens_cat_pts(1,2);
        
    case 'hslc'
        % Use HSLC as reference, fit other models to it
        true_params = params;
        true_params.sensory_info = sens_cat_pts(end,1);
        true_params.var_s = Model.getEvidenceVariance(sens_cat_pts(end,1));
        true_params.category_info = sens_cat_pts(end,2);
        true_params.p_match = sens_cat_pts(end,2);
end

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

%% Fit other models to the reference model

nModels = 4;
model_names = {'IS', 'gamma', 'ITB', 'ideal'};
% 'model_args' are arguments passed to @fit_helper below
model_args(1,:) = {'is',    {'log_temperature', 'log_lapse', 'prior_C', 'gamma', 'samples'}, true};
model_args(2,:) = {'ITB',   {'log_temperature', 'log_lapse', 'prior_C', 'gamma'},            true};
model_args(3,:) = {'ITB',   {'log_temperature', 'log_lapse', 'prior_C', 'gamma', 'bound'},   false};
model_args(4,:) = {'ideal', {'log_temperature', 'log_lapse', 'prior_C'},                     false};

for iModel=1:nModels
    fprintf('=== Fitting %s ===\n', model_names{iModel});
    bestfit{iModel} = fit_helper(true_params, model_args{iModel, :});
end

%% Model-comparison diagnostics

nReRuns = 100;
% ll(end,:) will be ll under true params
ll = zeros(nModels+1, nReRuns);
var_ll = zeros(nModels+1, nReRuns);
empty_prior = struct();

for iRun=nReRuns:-1:1
    fprintf('\tLL eval %d/%d\n', iRun, nReRuns);
    % Regenerate new random data
    true_params.seed = randi(1e9);
    data = Model.genDataWithParams(true_params);
    
    % Run true model and estimate its log likelihood on itself
    true_res = Model.runVectorized(true_params, data);
    [~, ll(end,iRun), var_ll(end,iRun)] = ...
        Fitting.choiceModelLogProb(true_params, empty_prior, data, true_res.choices);
    
    % Likewise for each of the other models
    for iModel=1:nModels
        [~, ll(iModel,iRun), var_ll(iModel,iRun)] = ...
            Fitting.choiceModelLogProb(bestfit{iModel}, empty_prior, data, true_res.choices);
    end
end

%% Compute LL _differences_ across reruns

% Each model saw the same data every time so their LLs are correlated across runs. The per-datum
% difference is more informative than the difference of overall performance.
ll_diff = ll - ll(end,:);
est_ll_diff = mean(ll_diff, 2);
sem_ll_diff = std(ll_diff, [], 2) ./ sqrt(nReRuns);

% For reference, compute the log likelihood under the null model
ll_null = -true_params.trials*log(2);

%% Plot result
figure; hold on;
bar(1:nModels, est_ll_diff(1:end-1));
errorbar(1:nModels, est_ll_diff(1:end-1), sem_ll_diff(1:end-1), 'ok');
% plot([1 5], ll_null-mean(ll(:,1)), '--k');
% text(1.5, ll_null-mean(ll(:,1))+5, 'random model');
set(gca, 'XTick', 1:nModels, 'XTickLabel', model_names);
grid on;
ylabel('\Delta LL from true model');

end

%% Helper function to get fit

function bestfit_params = fit_helper(true_params, model, fields, allow_gamma_neg)
uid = Model.getModelStringID(true_params);
fields_uid = strjoin(fields, '-');

for iRep=10:-1:1
    memo_file = fullfile(true_params.save_dir, ['bestfit-' model '-' uid '-' fields_uid '-' num2str(allow_gamma_neg) '-rep' num2str(iRep) '.mat']);
    [fit_params(iRep), fit_vals(iRep, :), nll(iRep), exitflag(iRep)] = LoadOrRun(@Model.fitModelToModel, ...
        {true_params, model, fields, allow_gamma_neg}, memo_file);
end

nll(exitflag <= 0) = inf;
[~, bestidx] = min(nll);
bestfit_params = fit_params(bestidx);

% figure;
% for i=1:5
%     for j=1:i
%         subplot(5,5,(i-1)*5+j); hold on;
%         plot(fit_vals(:,j),fit_vals(:,i), 'xk');
%         plot(Fitting.getParamsFields(true_params, fields{j}), Fitting.getParamsFields(true_params, fields{i}), 'or');
%
%         if j==1
%             ylabel(fields{i});
%         end
%         if i==5
%             xlabel(fields{j});
%         end
%     end
% end
end