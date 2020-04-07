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

% 'iModel' indexes the true model. 'jModel' indexes the fit model.
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
    both_data = {lshc_data{iModel}, hslc_data{iModel}};
    both_res = {lshc_res{iModel}, hslc_res{iModel}};
    
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

for ii=1:3*length(ground_truths)^2
    [iModel, jModel, iCondition] = ind2sub([length(ground_truths) length(ground_truths) 3], ii);
    base_params = ground_truths{jModel, 1};
    fields = ground_truths{jModel, 2};
    distribs = Fitting.defaultDistributions(fields, false, false);
    
    if jModel==1, continue; end
    
    switch iCondition
        case 1
            % LSHC fit
            uid = ['mhmap-' Model.getModelStringID(lshc_params{iModel}, true) '-lshc-' model_names{jModel} '-[' strjoin(fields, '-') ']'];
            [samples, lshc_fit_result{ii}, ~, sample_scores, metadata] = LoadOrRun(@Fitting.fitModelMH, ...
                {base_params, lshc_data{iModel}, lshc_res{iModel}.choices, distribs, struct('true_params', lshc_params{iModel})}, ...
                fullfile('../Precomputed/', [uid '.mat']));
            fprintf('FINISHED FIT %s TO %s [LSHC]\n', model_names{jModel}, model_names{iModel});
        case 2
            % HSLC fit
            uid = ['mhmap-' Model.getModelStringID(hslc_params{iModel}, true) '-hslc-' model_names{jModel} '-[' strjoin(fields, '-') ']'];
            [samples, hslc_fit_result{ii}, ~, sample_scores, metadata] = LoadOrRun(@Fitting.fitModelMH, ...
                {base_params, hslc_data{iModel}, hslc_res{iModel}.choices, distribs, struct('true_params', hslc_params{iModel})}, ...
                fullfile('../Precomputed/', [uid '.mat']));
            fprintf('FINISHED FIT %s TO %s [HSLC]\n', model_names{jModel}, model_names{iModel});
        case 3
            % Both-conditions fit
            both_data = [lshc_data{iModel}; hslc_data{iModel}];
            both_choices = [lshc_res{iModel}.choices; hslc_res{iModel}.choices];
            uid = ['mhmap-' Model.getModelStringID(hslc_params{iModel}, true) '-both-' model_names{jModel} '-[' strjoin(fields, '-') ']'];
            [samples, both_fit_result{ii}, ~, sample_scores, metadata] = LoadOrRun(@Fitting.fitModelMH, ...
                {base_params, both_data, both_choices, distribs, struct('true_params', [lshc_params{iModel} hslc_params{iModel}])}, ...
                fullfile('../Precomputed/', [uid '.mat']));
            fprintf('FINISHED FIT %s TO %s [BOTH]\n', model_names{jModel}, model_names{iModel});
    end
    
    clf;
    iplot = sample_scores.loglike > -sum([metadata.true_params.trials])*log(2);
    fprintf('%d of %d samples above null likelihiood threshold\n', sum(iplot), length(iplot));
    sz = 1./(1+exp(-zscore(sample_scores.loglike(iplot))));
    bestval = Fitting.getParamsFields(lshc_fit_result{ii}, fields);
    for iF=1:length(fields)
        for jF=iF:length(fields)
            subplot(length(fields),length(fields),(jF-1)*length(fields)+iF);
            hold on;
            scatter(samples(iplot,iF), samples(iplot,jF), 20, sample_scores.loglike(iplot), 'filled');
            plot(bestval(iF)*[1 1], [min(samples(:,jF)) max(samples(:,jF))], '--r');
            plot([min(samples(:,iF)) max(samples(:,iF))], bestval(jF)*[1 1], '--r');
            if iModel==jModel
                trueval = Fitting.getParamsFields(metadata.true_params(1), fields);
                plot(trueval(iF)*[1 1], [min(samples(:,jF)) max(samples(:,jF))], '--g');
                plot([min(samples(:,iF)) max(samples(:,iF))], trueval(jF)*[1 1], '--g');
            end
            
            xlabel(fields{iF}); ylabel(fields{jF});
        end
    end
    pause
end
