load(['../dscData/syntheticData_priorC5.mat']);
%%
nInner = 10;
fields = {'prior_C', 'gamma', 'samples', 'lapse', 'sensor_noise'};
LB = [0 0 1 0 0];
UB = [1 1 100 1 5];
PLB = LB;
PUB = UB;
x0 = [.5 .1 5 .01 1];
vbmc_options = vbmc('defaults');
vbmc_options.UncertaintyHandling = 'yes';
init_model_params = Model.newModelParams('model', 'is', ...
  'prior_C', 0.5, ...
  'gamma', 0.1, ...
  'samples', 5, ...
  'updates', 2, ...
  'lapse', 0.01);
 
% init_model_params = Fitting.sanitize(init_model_params);

likefn_args = {init_model_params , signals, sim_results.choices, nInner};
[VP, ~, ~, extflag] = LoadOrRun(@vbmc, ...
    [{@Fitting.choiceModelLogLikelihood2, x0, LB, UB, PLB, PUB, vbmc_options} likefn_args],...
    'vbmcRes');
%%
 