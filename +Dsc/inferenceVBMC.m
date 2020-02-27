clear all
clc
close all
%addpath(genpath('/Users/liushizhao/projects/APT/CB/confirmation-bias'));
%%
resultsFolder =  '../dscData/syntheticInferCB/';
i1 = 2;i2 = 2;i3 = 2;i4 = 2;i5 = 2;i6 =2;
dataName = sprintf('Trial%dpriorC%dgamma%dlapse%dsamples%dpmatch%dvars%d',...
    N_trial,i1,i2,i3,i4,i5,i6);
resultsName = dataName;
load([dataFolder,dataName])

% fields to infer
fields = {'prior_C', 'gamma', 'samples', 'lapse', 'sensor_noise'};
% upper and lower bound of each parameter
% is this the prior? is prior uniform distribution?
LB = [0 0 1 0 0];
UB = [1 1 100 1 5];
PLB = LB;
PUB = UB;
% initial value of parameters
x0 = [.5 .1 5 .01 1];


nInner = 10;
vbmc_options = vbmc('defaults');
vbmc_options.UncertaintyHandling = 'yes';
init_model_params = Model.newModelParams('model', 'is');
init_model_params = Fitting.sanitize(init_model_params);
likefn_args = {init_model_params , signals, sim_results.choices,nInner};


% define likelihood function
fun = @(x) llfun(x,fields,init_model_params,signals, sim_results.choices, nInner) ;
% vbmc inference
[VP, ~, ~, extflag] = vbmc(fun,x0,LB,UB,PLB,PUB,vbmc_options);
save([resultsFolder,resultsName],'VP','extflag')

function choiceLL = llfun(xvals,fields, model_params,signals, choices, nInner)
for iField=1:length(fields)
    % update inferred parameters for each iteration
    switch lower(fields{iField})
        case 'prior_c'
            model_params.prior_C = xvals(iField);
        case 'gamma'
            model_params.gamma = xvals(iField);
        case 'lapse'
            model_params.lapse = xvals(iField);
        case 'samples'
            model_params.samples = round(xvals(iField));
        case 'p_match'
            model_params.p_matcgh = xvals(iField);
        case 'var_s'
            model_params.var_s = xval(iField);
        case 'sensor_noise'
            sensor_noise = xvals(iField);
        case 'var_x'
            model_params.var_x = xvals(iField);
        case 'updates'
            model_params.updates = round(xvals(iField));
        otherwise
            error('Unrecognized / not implemented: %s', fields{iField});
    end
end
choiceLL = Fitting.choiceModelLogLikelihood(model_params, signals, choices, nInner);
end
