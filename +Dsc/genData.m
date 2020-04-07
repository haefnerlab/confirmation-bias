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
importantFields = {'prior_C', 'gamma', 'lapse', 'samples'};
defaultParams = Model.newModelParams('model', 'is');
dataFolder = '../dscData/syntheticDataCB/';
mkdir(dataFolder);

genDataInfo.N_trials_list = [800];
genDataInfo.prior_C_list = [0.3, 0.5, 0.7];
genDataInfo.gamma_list =  [0, 0.1, 0.3];
genDataInfo.lapse_list = [0, .1, .2];
genDataInfo.samples_list = [1,5,20];


save([dataFolder,'genDataInfo'],'genDataInfo','importantFields','defaultParams');
i0=1;
for i1 = 1:numel(genDataInfo.prior_C_list)
    for i2 = 1:numel(genDataInfo.gamma_list)
        for i3 = 1:numel(genDataInfo.lapse_list)
            for i4 = 1:numel(genDataInfo.samples_list)

                        params = Model.newModelParams('model', 'is', ...
                            'prior_C', genDataInfo.prior_C_list(i1),'gamma',genDataInfo.gamma_list(i2),...
                            'lapse',genDataInfo.lapse_list(i3),'samples',genDataInfo.samples_list(i4),...
                            'trials',genDataInfo.N_trials_list(i0));
                        [signals, categories,seed] = Model.genDataWithParams(params);
                        sim_results = Model.runVectorized(params, signals);
                        dataName = sprintf('Trial%dpriorC%dgamma%dlapse%dsamples%d',...
                            i0,i1,i2,i3,i4);
                        save([dataFolder,dataName],'signals','sim_results','params');

            end
        end
    end
end
