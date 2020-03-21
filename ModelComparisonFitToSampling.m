clear;

ps = 0.51:0.02:0.99;
true_params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'noise', 0, 'temperature', 0.1, 'trials', 1600, 'updates', 5, 'samples', 5);
true_params.save_dir = fullfile(pwd, 'tmp');
[~, sens_cat_pts] = Model.getThresholdPoints(ps, true_params, .7, 5);
rng('shuffle');
true_params.seed = randi(1e9);

lshc_params = true_params;
lshc_params.sensory_info = sens_cat_pts(1,1);
lshc_params.var_s = Model.getEvidenceVariance(sens_cat_pts(1,1));
lshc_params.category_info = sens_cat_pts(1,2);
lshc_params.p_match = sens_cat_pts(1,2);

hslc_params = true_params;
hslc_params.sensory_info = sens_cat_pts(end,1);
hslc_params.var_s = Model.getEvidenceVariance(sens_cat_pts(end,1));
hslc_params.category_info = sens_cat_pts(end,2);
hslc_params.p_match = sens_cat_pts(end,2);

%% Get best-fit ideal integrator.. just temperature, lapse, and prior fit
ideal_fit = fit_helper(lshc_params, 'ideal', {'log_temperature', 'log_lapse', 'prior_C'}, false);

%% Get best-fit bounded integrator
itb_fit = fit_helper(lshc_params, 'itb', {'log_temperature', 'log_lapse', 'prior_C', 'bound'}, false);

%% Get best-fit ITB-like integrator with CB term
cb_fit = fit_helper(lshc_params, 'itb', {'log_temperature', 'log_lapse', 'prior_C', 'gamma'}, true);

%% Get best-fit sampling model to itself
is_fit = fit_helper(lshc_params, 'is', {'log_temperature', 'log_lapse', 'prior_C', 'gamma', 'samples'}, true);

%% Model-comparison diagnostics

nReRuns = 100;
ll = zeros(5, nReRuns);
var_ll = zeros(5, nReRuns);
empty_prior = struct();

for iRun=nReRuns:-1:1
    disp(iRun);
    % Regenerate new random data
    lshc_params.seed = randi(1e9);
    data = Model.genDataWithParams(lshc_params);
    
    % Run true model and estimate its log likelihood on itself
    true_res = Model.runVectorized(lshc_params, data);
    [~, ll(1,iRun), var_ll(1,iRun)] = Fitting.choiceModelLogProb(lshc_params, empty_prior, data, true_res.choices);
    % Likewise for model that was fit from the 'correct' model family
    [~, ll(2,iRun), var_ll(2,iRun)] = Fitting.choiceModelLogProb(is_fit, empty_prior, data, true_res.choices);
    % Likewise for CB-like model
    [~, ll(3,iRun), var_ll(3,iRun)] = Fitting.choiceModelLogProb(cb_fit, empty_prior, data, true_res.choices);
    % Estimate log likelihood under ideal model
    [~, ll(4,iRun), var_ll(4,iRun)] = Fitting.choiceModelLogProb(ideal_fit, empty_prior, data, true_res.choices);
    % Likewise for ITB
    [~, ll(5,iRun), var_ll(5,iRun)] = Fitting.choiceModelLogProb(itb_fit, empty_prior, data, true_res.choices);
end

est_ll = mean(ll, 2);
est_ll_sem = mean(var_ll, 2) / sqrt(nReRuns);

% Compute difference in LL per run; each model saw the same data every time so their LLs are
% correlated across runs. The per-datum difference is more informative than the difference of
% overall performance.
ll_diff = ll - ll(1, :);
est_ll_diff = mean(ll_diff, 2);
sem_ll_diff = std(ll_diff, [], 2) ./ sqrt(nReRuns);

% For reference, compute the log likelihood under the null model
ll_null = -true_params.trials*log(2);

%% Plot result
figure; hold on;
bar(1:5, est_ll_diff);
errorbar(1:5, est_ll_diff, sem_ll_diff, 'ok');
% plot([1 5], ll_null-mean(ll(:,1)), '--k');
% text(1.5, ll_null-mean(ll(:,1))+5, 'random model');
set(gca, 'XTick', 1:5, 'XTickLabel', {'true model', 'IS fit', 'CB', 'ideal', 'ITB'});
ylabel('log likelihood');

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