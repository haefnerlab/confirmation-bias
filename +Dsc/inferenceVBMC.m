clear
clc
close all
addpath(genpath('../dscData/'));
%%
i = 5;
load(['syntheticData_priorC',num2str(i)]);

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
    'gamma', 0.1, ...
    'samples', 5, ...
    'updates', 2, ...
    'lapse', 0.01);
init_model_params = Fitting.sanitize(init_model_params);
likefn_args = {init_model_params , signals, sim_results.choices,nInner};
% [VP, ~, ~, extflag] = LoadOrRun(@vbmc, ...
%         [{@Fitting.choiceModelLogLikelihood, x0, LB, UB, PLB, PUB, vbmc_options},likefn_args],...
%         'vbmcRes');
%%
% llfun = @Fitting.choiceModelLogLikelihood;
fun = @(x) llfun(x,fields,init_model_params,signals, sim_results.choices, nInner) ;
% Fitting.choiceModelLogLikelihood(init_model_params , signals, sim_results.choices, nInner)
[VP, ~, ~, extflag] = vbmc(fun,x0,LB,UB,PLB,PUB,vbmc_options);

Xsamp = vbmc_rnd(VP, 1e5);
[fig, ax] = cornerplot(Xsamp, fields, [], [LB; UB]);
statuses = {'not converged', 'converged'};

disp(statuses{extflag+1});



function choicell = llfun(xvals,fields, model_params, signals, choices, nInner)
    for iField=1:length(fields)
        switch lower(fields{iField})
            case 'prior_c'
                model_params.prior_C = xvals(iField);
            case 'lapse'
                model_params.lapse = xvals(iField);
            case 'gamma'
                model_params.gamma = xvals(iField);
            case 'sensor_noise'
                model_params.sensor_noise = xvals(iField);
            case 'var_x'
                model_params.var_x = xvals(iField);
            case 'updates'
                model_params.updates = round(xvals(iField));
            case 'samples'
                model_params.samples = round(xvals(iField));
            otherwise
                error('Unrecognized / not implemented: %s', fields{iField});
        end
    end
    choicell = Fitting.choiceModelLogLikelihood(model_params, signals, choices, nInner);
end