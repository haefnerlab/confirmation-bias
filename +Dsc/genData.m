clear all
clc
close all
% addpath(genpath('/Users/liushizhao/projects/APT/CB/confirmation-bias'));
%%
% generate data with chosen parameters
fields_domains = {'var_s',  logspace(-2, 1), true;
    'var_x',  logspace(-2, 1), true;
    'prior_C',  linspace(0, 1), false;
    'p_match',  linspace(0, 1), false;
    'gamma',  linspace(0, 1), false;
    'noise', linspace(0, 2), false;
    'updates',  1:100, false;
    'lapse', linspace(.001, 1), false};
important_fields = {'prior_C', 'gamma', 'lapse', 'samples', 'p_match', 'var_s'};
prior_C_list = linspace(0, 1,10);
for i = 1:numel(prior_C_list)
    params = Model.newModelParams('model', 'is', ...
        'prior_C', prior_C_list(i),'trials',10);
    [signals, categories,seed] = Model.genDataWithParams(params);
    sim_results = Model.runVectorized(params, signals);

    save(['../dscData/syntheticData_priorC',num2str(i)],'signals','sim_results','params');
end
