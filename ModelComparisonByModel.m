ground_truths = {
    Model.newModelParams('model', 'is', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1, 'gamma', 0.1, 'samples', 5, 'updates', 5, 'seed', 18863457)
    Model.newModelParams('model', 'itb', 'trials', 1000, 'temperature', 0.1, 'bound', 1.2, 'gammafun', @(ci,si) 1-ci, 'noise', 0.35, 'updates', 1, 'seed', 18863457);
};


lshc_params = cell(size(ground_truths));
lshc_data = cell(size(ground_truths));
lshc_res = cell(size(ground_truths));

hslc_params = cell(size(ground_truths));
hslc_data = cell(size(ground_truths));
hslc_res = cell(size(ground_truths));

both_params = cell(size(ground_truths));
both_data = cell(size(ground_truths));
both_res = cell(size(ground_truths));

for iModel=1:length(ground_truths)
    true_params = ground_truths{iModel};
    sens_cat_pts = Model.getThresholdPoints(0.51:0.02:0.99, true_params, .7, 5);
    
    % LSHC ground truth
    lshc_params{iModel} = Model.setCategorySensoryInfo(true_params, sens_cat_pts(1,2), sens_cat_pts(1,1));
    lshc_data{iModel} = Model.genDataWithParams(lshc_params{iModel});
    lshc_res{iModel} = Model.runVectorized(lshc_params{iModel}, lshc_data{iModel});
    
    % HSLC ground truth
    hslc_params{iModel} = Model.setCategorySensoryInfo(true_params, sens_cat_pts(end,2), sens_cat_pts(end,1));
    hslc_data{iModel} = Model.genDataWithParams(hslc_params{iModel});
    hslc_res{iModel} = Model.runVectorized(hslc_params{iModel}, hslc_data{iModel});
    
    % Combine 'em
    both_params{iModel} = [lshc_params{iModel}, hslc_params{iModel}];
    both_data{iModel} = {lshc_data{iModel}, hslc_data{iModel}};
    both_res{iModel} = [lshc_res{iModel}, hslc_res{iModel}];
    
    % DEBUGGING: plot PKs
    subplot(2,length(ground_truths),iModel);
    [pk, ~, err] = CustomRegression.PsychophysicalKernel(lshc_data{iModel}, lshc_res{iModel}.choices == +1, 1, 0, 100, 1);
    errorbar(pk, err);
    title([true_params.model ' [LSHC]']);
    subplot(2,length(ground_truths),iModel+length(ground_truths));
    [pk, ~, err] = CustomRegression.PsychophysicalKernel(hslc_data{iModel}, hslc_res{iModel}.choices == +1, 1, 0, 100, 1);
    errorbar(pk, err);
    title([true_params.model ' [HSLC]']);
end

%%

parfor ii=1:3*length(ground_truths)
    [iModel, iCondition] = ind2sub([length(ground_truths) 3], ii);
        
    switch iCondition
        case 1
            prefix = ['gt-' Model.getModelStringID(lshc_params{iModel}, true) '-lshc'];
            [aic, model_info, sampleses] = ModelComparison(lshc_params{iModel}, lshc_data{iModel}, lshc_res{iModel}.choices, prefix);
        case 2
            prefix = ['gt-' Model.getModelStringID(hslc_params{iModel}, true) '-hslc'];
            [aic, model_info, sampleses] = ModelComparison(hslc_params{iModel}, hslc_data{iModel}, hslc_res{iModel}.choices, prefix);
        case 3
            prefix = ['gt-' Model.getModelStringID(lshc_params{iModel}, true) '-both'];
            [aic, model_info, sampleses] = ModelComparison(both_params{iModel}, both_data{iModel}, {both_res{iModel}.choices}, prefix);
    end
end
