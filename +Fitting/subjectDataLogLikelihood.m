function [loglike, model_params, data_translated] = subjectDataLogLikelihood(xvals, fields, model_params, SubjectData, nInner)

% Set defaults
sensor_noise = 0.1;

for iField=1:length(fields)
    switch lower(fields{iField})
        case 'prior_c'
            model_params.prior_C = xvals(iField);
        case 'lapse'
            model_params.lapse = xvals(iField);
        case 'gamma'
            model_params.gamma = xvals(iField);
        case 'sensor_noise'
            sensor_noise = xvals(iField);
        case 'var_x'
            model_params.var_x = xvals(iField);
        case 'updates'
            model_params.updates = round(xvals(iField));
        case 'samples'
            model_params.samples = round(xvals(iField));
        otherwise
            error('Unrecognized / not implemented: %s', fields{iField});
    end
end

model_params = Fitting.sanitize(model_params);

unsigned_ratio = max(SubjectData.true_ratio, 1-SubjectData.true_ratio);
[kvals, kcounts] = unique_counts(SubjectData.noise);
[pvals, pcounts] = unique_counts(unsigned_ratio);

[~, sensory_llos, empirical_sensory_info] = bpg.signalToLLO(...
        SubjectData.ideal_frame_signals, [kvals kcounts/sum(kcounts)], [pvals pcounts/sum(pcounts)], sensor_noise);

% Assume optimal (average) value of var_s and p_match (TODO - make these parameters)
opt_sensory_info = mean(empirical_sensory_info(:));
opt_category_info = dot(pvals, pcounts/sum(pcounts));
model_params.var_s = Model.getEvidenceVariance(opt_sensory_info);
model_params.p_match = opt_category_info;

% Translate from observed LLO to model inputs
data_translated = Model.lloToEvidence(model_params.var_s, sensory_llos);

% Compute 'choice' log likelihood of model on pseudo-data as fed into the model
loglike = Fitting.choiceModelLogLikelihood(model_params, data_translated, SubjectData.choice == +1, nInner);

end