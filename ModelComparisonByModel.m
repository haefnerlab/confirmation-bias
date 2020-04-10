ground_truths = {
    Model.newModelParams('model', 'is', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1, 'gamma', 0.1, 'samples', 5, 'updates', 5, 'seed', 18863457)
    Model.newModelParams('model', 'itb', 'trials', 1000, 'temperature', 0.1, 'bound', 1.2, 'gammafun', @(ci,si) 1-ci, 'noise', 0.35, 'updates', 1, 'seed', 18863457);
};

% Names copied from @ModelComparison
model_names = {'is', 'vb', 'itb', 'itb-gamma', 'ideal'};


lshc_params = cell(size(ground_truths));
lshc_data = cell(size(ground_truths));
lshc_res = cell(size(ground_truths));

hslc_params = cell(size(ground_truths));
hslc_data = cell(size(ground_truths));
hslc_res = cell(size(ground_truths));

both_params = cell(size(ground_truths));
both_data = cell(size(ground_truths));
both_res = cell(size(ground_truths));

for iTruth=1:length(ground_truths)
    true_params = ground_truths{iTruth};
    sens_cat_pts = Model.getThresholdPoints(0.51:0.02:0.99, true_params, .7, 5);
    
    % LSHC ground truth
    lshc_params{iTruth} = Model.setCategorySensoryInfo(true_params, sens_cat_pts(1,2), sens_cat_pts(1,1));
    lshc_data{iTruth} = Model.genDataWithParams(lshc_params{iTruth});
    lshc_res{iTruth} = Model.runVectorized(lshc_params{iTruth}, lshc_data{iTruth});
    
    % HSLC ground truth
    hslc_params{iTruth} = Model.setCategorySensoryInfo(true_params, sens_cat_pts(end,2), sens_cat_pts(end,1));
    hslc_data{iTruth} = Model.genDataWithParams(hslc_params{iTruth});
    hslc_res{iTruth} = Model.runVectorized(hslc_params{iTruth}, hslc_data{iTruth});
    
    % Combine 'em
    both_params{iTruth} = [lshc_params{iTruth}, hslc_params{iTruth}];
    both_data{iTruth} = {lshc_data{iTruth}, hslc_data{iTruth}};
    both_res{iTruth} = [lshc_res{iTruth}, hslc_res{iTruth}];
    
    % DEBUGGING: plot PKs
    subplot(2,length(ground_truths),iTruth);
    [pk, ~, err] = CustomRegression.PsychophysicalKernel(lshc_data{iTruth}, lshc_res{iTruth}.choices == +1, 1, 0, 100, 1);
    errorbar(pk, err);
    title([true_params.model ' [LSHC]']);
    subplot(2,length(ground_truths),iTruth+length(ground_truths));
    [pk, ~, err] = CustomRegression.PsychophysicalKernel(hslc_data{iTruth}, hslc_res{iTruth}.choices == +1, 1, 0, 100, 1);
    errorbar(pk, err);
    title([true_params.model ' [HSLC]']);
end

%%

parfor ii=1:3*length(ground_truths)*length(model_names)
    [iTruth, iModel, iCondition] = ind2sub([length(ground_truths) length(model_names) 3], ii);

    switch iCondition
        case 1
            prefix = ['gt-' Model.getModelStringID(lshc_params{iTruth}, true) '-lshc'];
            [aic, ~, model_info] = ModelComparison(lshc_params{iTruth}, lshc_data{iTruth}, lshc_res{iTruth}.choices, false, prefix, model_names(iModel));
        case 2
            prefix = ['gt-' Model.getModelStringID(hslc_params{iTruth}, true) '-hslc'];
            [aic, ~, model_info] = ModelComparison(hslc_params{iTruth}, hslc_data{iTruth}, hslc_res{iTruth}.choices, false, prefix, model_names(iModel));
        case 3
            prefix = ['gt-' Model.getModelStringID(lshc_params{iTruth}, true) '-both'];
            [aic, ~, model_info] = ModelComparison(both_params{iTruth}, both_data{iTruth}, {both_res{iTruth}.choices}, false, prefix, model_names(iModel));
    end
end
