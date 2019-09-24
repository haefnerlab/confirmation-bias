function log_likelihood = pkModelLogLikelihood(params, pk_mean, pk_var)

% Clear seed for random simulation each time
params.seed = [];
params.var_s = params.var_s_per_sample * params.samples;
data = Model.genDataWithParams(params);
results = Model.runVectorized(params, data);
[model_pk, ~, est_std] = CustomRegression.PsychophysicalKernel(data, results.choices == +1, 0, 0, 0, 0);
net_var = pk_var(:) + est_std(:).^2;
pk_diffs = pk_mean(:) - model_pk(:);
log_likelihood = -1/2 * sum(pk_diffs.^2 ./ net_var);

end