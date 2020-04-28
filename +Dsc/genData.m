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

                        base_params = Model.newModelParams('model', 'is', ...
                            'prior_C', genDataInfo.prior_C_list(i1),'gamma',genDataInfo.gamma_list(i2),...
                            'lapse',genDataInfo.lapse_list(i3),'samples',genDataInfo.samples_list(i4),...
                            'trials',genDataInfo.N_trials_list(i0));
                        % Find points along the 70% performance curve, ordered from lowest to highest sensory info
                        performance_level = 0.7;
                        si_ci_grid = 0.51:0.02:0.99;
                        sens_cat_pts = Model.getThresholdPoints(si_ci_grid, base_params, performance_level, 5);
                        % Set LSHC params to be equal to lowest of the above SI values at 70% performance
                        lshc_params = Model.setCategorySensoryInfo(base_params, sens_cat_pts(1,2), sens_cat_pts(1,1));
                        % Likewise for HSLC using highest of the SI values above
                        hslc_params = Model.setCategorySensoryInfo(base_params, sens_cat_pts(end,2), sens_cat_pts(end,1));
                        
                        lshc_signals = Model.genDataWithParams(lshc_params);
                        hslc_signals = Model.genDataWithParams(hslc_params);
                        
                        lshc_sim_results = Model.runVectorized(lshc_params, lshc_signals);
                        hslc_sim_results = Model.runVectorized(hslc_params, hslc_signals);
                        dataName = sprintf('Trial%dpriorC%dgamma%dlapse%dsamples%d',...
                            i0,i1,i2,i3,i4);
                        save([dataFolder,dataName],'hslc_signals','lshc_signals','lshc_sim_results',...
                            'hslc_sim_results','lshc_params','hslc_params');

            end
        end
    end
end
