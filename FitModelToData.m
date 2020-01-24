% function [bestfit_params, loglike] = FitModelToData(SubjectData, varargin)

init_model_params = Model.newModelParams('model', 'is', ...
'gamma', 0.1, ...
'samples', 5, ...
'updates', 2, ...
'lapse', 0.01);

nInner = 10;

fields = {'prior_C', 'gamma', 'samples', 'lapse', 'sensor_noise'};
LB = [0 0 1 0 0];
UB = [1 1 100 1 5];
PLB = LB;
PUB = UB;

x0 = [.5 .1 5 .01 1];

likefn_args = {fields, init_model_params, ThresholdData, nInner};

vbmc_options = vbmc('defaults');
vbmc_options.UncertaintyHandling = 'yes';
[VP, ELBO, ELBO_SD, EXITFLAG] = vbmc(@Fitting.subjectDataLogLikelihood, ...
    x0, LB, UB, PLB, PUB, vbmc_options, likefn_args{:});


% end