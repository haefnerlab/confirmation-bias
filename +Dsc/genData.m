clear
clc
close all
% addpath(genpath('/Users/liushizhao/projects/APT/CB/confirmation-bias'));
%%
% generate data with chosen parameters
% fields_domains = {'var_s',  logspace(-2, 1), true;
%     'var_x',  logspace(-2, 1), true;
%     'prior_C',  linspace(0, 1), false;
%     'p_match',  linspace(0, 1), false;
%     'gamma',  linspace(0, 1), false;
%     'noise', linspace(0, 2), false;
%     'updates',  1:100, false;
%     'lapse', linspace(.001, 1), false};
importantFields = {'prior_C', 'gamma', 'lapse', 'samples', 'p_match', 'var_s'};
defaultParams = Model.newModelParams('model', 'is');
dataFolder = '../dscData/syntheticDataCB/';
mkdir(dataFolder);

genDataInfo.N_trials_list = [10, 200, 800];
genDataInfo.prior_C_list = [0.1, 0.5, 0.9];
genDataInfo.gamma_list =  [0, 0.2, 0.4, 0.6];
genDataInfo.lapse_list = [0, .1, .2];
genDataInfo.samples_list = [1,5,10,50];
genDataInfo.p_match_list = [0, 0.5, 1];
genDataInfo.var_s_list = logspace(-2,1,3);

save([dataFolder,'genDataInfo'],'genDataInfo','importantFields','defaultParams');
for i0 = 1:numel(genDataInfo.N_trials_list)
for i1 = 1:numel(genDataInfo.prior_C_list)
    for i2 = 1:numel(genDataInfo.gamma_list)
        for i3 = 1:numel(genDataInfo.lapse_list)
            for i4 = 1:numel(genDataInfo.samples_list)
                for i5 = 1:numel(genDataInfo.p_match_list)
                    for i6 = 1:numel(genDataInfo.var_s_list)
                        params = Model.newModelParams('model', 'is', ...
                            'prior_C', genDataInfo.prior_C_list(i1),'gamma',genDataInfo.gamma_list(i2),...
                            'lapse',genDataInfo.lapse_list(i3),'samples',genDataInfo.samples_list(i4),...
                            'p_match',genDataInfo.p_match_list(i5),'var_s',genDataInfo.var_s_list(i6),...
                            'trials',genDataInfo.N_trials_list(i0));
                        [signals, categories,seed] = Model.genDataWithParams(params);
                        sim_results = Model.runVectorized(params, signals);
                        dataName = sprintf('Trial%dpriorC%dgamma%dlapse%dsamples%dpmatch%dvars%d',...
                            i0,i1,i2,i3,i4,i5,i6);
                        save([dataFolder,dataName],'signals','sim_results','params');
                    end
                end
            end
        end
    end
end
end