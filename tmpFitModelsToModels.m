model_params_fields = {
    Model.newModelParams('model', 'is', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1, 'gamma', 0.1, 'samples', 5, 'updates', 5), {'prior_C', 'temperature', 'lapse', 'gamma', 'samples'}
    Model.newModelParams('model', 'vb-czx', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1, 'gamma', 0.1, 'updates', 5, 'step_size', 0.05), {'prior_C', 'temperature', 'lapse', 'gamma', 'step_size'}
    Model.newModelParams('model', 'itb', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1, 'gamma', 0.1, 'updates', 1, 'bound', 1), {'prior_C', 'temperature', 'lapse', 'gamma', 'bound'}
    Model.newModelParams('model', 'ideal', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1), {'prior_C', 'temperature', 'lapse'}
    };

model_names = {'is', 'vb-czx', 'itb', 'ideal'};

rng(2745393, 'twister');
seeds = randi(1e9, length(model_names), 1);

lshc_params = cell(4,1);
lshc_res = cell(4,1);
lshc_fit_result = cell(4,4);
hslc_params = cell(4,1);
hslc_res = cell(4,1);
hslc_fit_result = cell(4,4);

% 'iModel' indexes the true model. 'jModel' indexes the fit model.
for iModel=1:length(model_params_fields)
    true_params = model_params_fields{iModel,1};
    true_params.save_dir = 'tmp';
    sens_cat_pts = Model.getThresholdPoints(0.51:0.02:0.99, true_params, .7, 5);
    
    % LSHC ground truth
    lshc_params{iModel} = true_params;
    lshc_params{iModel}.sensory_info = sens_cat_pts(1,1);
    lshc_params{iModel}.var_s = Model.getEvidenceVariance(sens_cat_pts(1,1));
    lshc_params{iModel}.category_info = sens_cat_pts(1,2);
    lshc_params{iModel}.p_match = sens_cat_pts(1,2);
    
    lshc_params{iModel}.seed = seeds(iModel);
    lshc_data = Model.genDataWithParams(lshc_params{iModel});
    lshc_res{iModel} = Model.runVectorized(lshc_params{iModel}, lshc_data);
    
    % HSLC ground truth
    hslc_params{iModel} = true_params;
    hslc_params{iModel}.sensory_info = sens_cat_pts(1,1);
    hslc_params{iModel}.var_s = Model.getEvidenceVariance(sens_cat_pts(1,1));
    hslc_params{iModel}.category_info = sens_cat_pts(1,2);
    hslc_params{iModel}.p_match = sens_cat_pts(1,2);
    
    hslc_params{iModel}.seed = seeds(iModel);
    hslc_data = Model.genDataWithParams(hslc_params{iModel});
    hslc_res{iModel} = Model.runVectorized(hslc_params{iModel}, hslc_data);
end

parfor ii=1:length(model_params_fields)^2
    [iModel, jModel] = ind2sub([length(model_params_fields) length(model_params_fields)], ii);
    base_params = model_params_fields{jModel, 1};
    fields = model_params_fields{jModel, 2};
    distribs = Fitting.defaultDistributions(fields, false, false);
    
    % LSHC fit
    uid = ['mhmap-' Model.getModelStringID(lshc_params{iModel}, true) '-lshc-' model_names{jModel} '-[' strjoin(fields, '-') ']'];
    [lshc_fit_result{ii}, ~, ~, ~] = LoadOrRun(@Fitting.fitModelMHMAP, ...
        {base_params, lshc_data, lshc_res{iModel}.choices, distribs, struct('true_params', lshc_params{iModel})}, ...
        fullfile('../Precomputed/', [uid '.mat']));
    fprintf('FINISHED FIT %s TO %s [LSHC]\n', model_names{jModel}, model_names{iModel});
    
%     clf;
%     sz = 1./(1+exp(-zscore(scores)));
%     bestval = Fitting.getParamsFields(fit_result{iModel,jModel}, fields);
%     for iF=1:length(fields)
%         for jF=iF:length(fields)
%             subplot(length(fields),length(fields),(jF-1)*length(fields)+iF);
%             hold on;
%             scatter(samples(:,iF), samples(:,jF), sz*30, 'filled');
%             plot(bestval(iF)*[1 1], [min(samples(:,jF)) max(samples(:,jF))], '--r');
%             plot([min(samples(:,iF)) max(samples(:,iF))], bestval(jF)*[1 1], '--r');
%             if iModel==jModel
%                 trueval = Fitting.getParamsFields(metadata.true_params, fields);
%                 plot(trueval(iF)*[1 1], [min(samples(:,jF)) max(samples(:,jF))], '--g');
%                 plot([min(samples(:,iF)) max(samples(:,iF))], trueval(jF)*[1 1], '--g');
%             end
%             
%             xlabel(fields{iF}); ylabel(fields{jF});
%         end
%     end
    
    % HSLC fit
    uid = ['mhmap-' Model.getModelStringID(hslc_params{iModel}, true) '-hslc-' model_names{jModel} '-[' strjoin(fields, '-') ']'];
    [hslc_fit_result{ii}, ~, ~, ~] = LoadOrRun(@Fitting.fitModelMHMAP, ...
        {base_params, hslc_data, hslc_res{iModel}.choices, distribs, struct('true_params', hslc_params{iModel})}, ...
        fullfile('../Precomputed/', [uid '.mat']));
    fprintf('FINISHED FIT %s TO %s [HSLC]\n', model_names{jModel}, model_names{iModel});
end