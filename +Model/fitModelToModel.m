function [fit_params, fit_vals, nll, exitflag] = fitModelToModel(true_params, modeltofit, fields, allow_gamma_neg)
base_params = true_params;
base_params.model = lower(modeltofit);
base_params.allow_gamma_neg = allow_gamma_neg;
distribs = Fitting.defaultDistributions(fields, true, allow_gamma_neg);

data = Model.genDataWithParams(true_params);
true_res = Model.runVectorized(true_params, data);

%% Do fit
LB = cellfun(@(f) distribs.(f).lb, fields);
UB = cellfun(@(f) distribs.(f).ub, fields);
PLB = cellfun(@(f) distribs.(f).plb, fields);
PUB = cellfun(@(f) distribs.(f).pub, fields);

options = bads('defaults');
options.Display = 'iter';
options.UncertaintyHandling = true;
options.NonlinearScaling = false;

%% Round 0: evaluate draws from the prior and choose best draw as initialization

% Draw a batch of initial values from the prior and evaluate them
nInitDraw = 100;
for iField=length(fields):-1:1
    init_draw(:,iField) = distribs.(fields{iField}).priorrnd([nInitDraw 1]);
end

for iDraw=nInitDraw:-1:1
    [init_post(iDraw), ~, init_var(iDraw)] = Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, init_draw(iDraw, :)), ...
        distribs, data, true_res.choices);
end
init_upper_conf = init_post + 3*sqrt(init_var);

[~, best_init] = max(init_upper_conf);
init_vals = init_draw(best_init, :);

%% Round 1: run BADS from best initialization using appropriate noise size

% Estimate size of noise at that initial value and inflate it for stability of fit
[~,~,est_var] = Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, init_vals), ...
    distribs, data, true_res.choices);
options.NoiseSize = sqrt(2*est_var);

% Run optimization
[fit_vals, ~, ~] = bads(...
    @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, x), distribs, data, true_res.choices), ...
    init_vals, LB, UB, PLB, PUB, [], options);

%% Round 2: re-estimate noise and run again from the best point

[~,~,est_var] = Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, fit_vals), ...
    distribs, data, true_res.choices);
options.NoiseSize = sqrt(est_var);

% Re-run optimization starting from best point in round 1
[fit_vals, nll, exitflag] = bads(...
    @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, x), distribs, data, true_res.choices), ...
    fit_vals, LB, UB, PLB, PUB, [], options);

%% Return best fit
fit_params = Fitting.setParamsFields(base_params, fields, fit_vals);
end