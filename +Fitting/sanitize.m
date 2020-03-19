function params = sanitize(params)
params.prior_C = clip(params.prior_C, 0, 1);
params.lapse = clip(params.lapse, 1e-9, 1-1e-9);
params.noise = clip(params.noise, 1e-9, inf);
params.temperature = clip(params.temperature, 1e-9, inf);
if params.allow_gamma_neg
	params.gamma = clip(params.gamma, -1, 1);
else
	params.gamma = clip(params.gamma, 0, 1);
end
% During fitting of IS model, var_s and samples are highly correlated, so we reparameterize var_s as
% variance per sample.
if isfield(params, 'var_s_per_sample')
    params.var_s_per_sample = clip(params.var_s_per_sample, 1e-9, inf);
    params.var_s = params.var_s_per_sample * params.samples;
else
    params.var_s = clip(params.var_s, 1e-9, inf);
    params.var_s_per_sample = params.var_s / params.samples;
end
params.var_x = clip(params.var_x, 1e-9, inf);
params.updates = clip(round(params.updates), 1, inf);
params.samples = clip(round(params.samples), 1, inf);
params.bound = clip(params.bound, 1e-9, inf);

end

function val = clip(val, lo, hi)
val = min(max(val, lo), hi);
end