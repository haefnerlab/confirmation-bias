clear
clc
close all
%addpath(genpath('/Users/liushizhao/projects/APT/CB/confirmation-bias'));
%%
clear;
dataFolder = '../dscData/syntheticDataCB/';
resultsFolder =  '../dscData/syntheticInferCBresults/';
load([dataFolder,'genDataInfo', '.mat'], 'genDataInfo');
fieldsToInfer = {'prior_C', 'gamma', 'samples', 'lapse'};
% for i0 = 1:numel(genDataInfo.N_trials_list)
% for i1 = 1:numel(genDataInfo.prior_C_list)
%     for i2 = 1:numel(genDataInfo.gamma_list)
%         for i3 = 1:numel(genDataInfo.lapse_list)
%             for i4 = 1:numel(genDataInfo.samples_list)
%                 for i5 = 1:numel(genDataInfo.p_match_list)
%                     for i6 = 1:numel(genDataInfo.var_s_list)
%                         [VP, extflag, dataName] = runVBMCInferenceFunc(dataFolder, fieldsToInfer, i0, i1, i2, i3, i4, i5, i6);
%                         parsave([resultsFolder,dataName],VP,extflag,fieldsToInfer)
%                     end
%                 end
%             end
%         end
%     end
% end
% end
A = [1, 2, 3]';
ma = size(A, 1);
[a,b,c, d, e, f, g]=ndgrid(1:ma,1:ma,1:ma, 1:ma, 1:ma, 1:ma, 1:ma);
product = [A(a,:),A(b,:),A(c,:),A(d, :), A(e, :), A(f, :), A(g, :)];
parfor (i = 1:size(product, 1))
    [VP, extflag, dataName] = runVBMCInferenceFunc(dataFolder, fieldsToInfer, ...
                                product(i, 1), product(i, 2), product(i, 3), ...
                                product(i, 4), product(i, 5), product(i, 6), ...
                                product(i, 7));
    parsave([resultsFolder,dataName],VP,extflag,fieldsToInfer)
    if mod(i, 10) == 0
        try
            fid = fopen('output_file.txt', 'at+');
            fprintf(fid, "Done with: %d samples", i/10);
            fclose(fid);
            disp(["Done with, ", num2str(i/10), " samples"]);
        catch
            disp("That didn't work");
        end
    end
end
% [VP, extflag, dataName] = runVBMCInferenceFunc(dataFolder, fieldsToInfer, i0, i1, i2, i3, i4, i5, i6);
% Xsamp = vbmc_rnd(VP, 1e5);
% [fig, ax] = cornerplot(Xsamp, fields, [], [LB; UB]);
% save([resultsFolder,dataName],'VP','extflag','fieldsToInfer')
% save the figure here?
% save the inferred variables here
function parsave(fname, VP, extflag, fieldsToInfer)
  save(fname,'VP','extflag','fieldsToInfer');
end
function [VP, extflag, dataName] = runVBMCInferenceFunc(dataFolder, fieldsToInfer, ...
        i_trials, i_prior, i_gamma, i_lapse, i_samples, i_p_match, i_var_s)
    dataName = sprintf('Trial%dpriorC%dgamma%dlapse%dsamples%dpmatch%dvars%d',...
        i_trials,i_prior,i_gamma, i_lapse, i_samples, i_p_match, i_var_s);
    load([dataFolder,dataName, '.mat']);
    % fields to infer
    % upper and lower bound of each parameter
    % is this the prior? is prior uniform distribution?
    LB = [0 0 1 0 ];
    UB = [1 1 100 1];
    PLB = LB;
    PUB = UB;
    % initial value of parameters
    x0 = [.5 .1 5 .01];
    nInner = 10;
    vbmc_options = vbmc('defaults');
    vbmc_options.UncertaintyHandling = 'yes';
    vbmc_options.Display = 'final';
    init_model_params = Model.newModelParams('model', 'is');
    init_model_params = Fitting.sanitize(init_model_params);
    likefn_args = {init_model_params , signals, sim_results.choices,nInner};
    % define likelihood function
    fun = @(x) llfun(x,fieldsToInfer,init_model_params,signals, sim_results.choices, nInner) ;
    % vbmc inference
    [VP, ~, ~, extflag] = vbmc(fun,x0,LB,UB,PLB,PUB,vbmc_options);
end
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
%         case 'sensor_noise'
%             sensor_noise = xvals(iField);
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