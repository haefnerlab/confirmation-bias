function params = sanitize(params)

params.prior_C = clip(params.prior_C, 0, 1);
params.lapse = clip(params.lapse, 0, 1);
params.gamma = clip(params.gamma, 0, 1);
if isfield(params, 'var_s_per_sample')
    params.var_s_per_sample = clip(params.var_s_per_sample, 0, inf);
    params.var_s = params.var_s_per_sample * params.samples;
else
    params.var_s = clip(params.var_s, 0, inf);
    params.var_s_per_sample = params.var_s / params.samples;
end
params.var_x = clip(params.var_x, 0, inf);
params.updates = clip(round(params.updates), 1, inf);
params.samples = clip(round(params.samples), 1, inf);

end

function val = clip(val, lo, hi)
val = min(max(val, lo), hi);
end