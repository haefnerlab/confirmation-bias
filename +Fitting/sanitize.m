function params = sanitize(params)
fields = fieldnames(params);
for iF=1:length(fields)
    if startsWith(fields{iF}, 'prior_C', 'IgnoreCase', true)
        params.(fields{iF}) = clip(params.(fields{iF}), 0, 1);
    elseif startsWith(fields{iF}, 'lapse', 'IgnoreCase', true)
        params.(fields{iF}) = clip(params.(fields{iF}), 1e-9, 1-1e-9);
    elseif startsWith(fields{iF}, 'noise', 'IgnoreCase', true)
        params.(fields{iF}) = clip(params.(fields{iF}), 1e-9, inf);
    elseif startsWith(fields{iF}, 'temperature', 'IgnoreCase', true)
        params.(fields{iF}) = clip(params.(fields{iF}), 1e-9, inf);
    elseif startsWith(fields{iF}, 'gamma', 'IgnoreCase', true)
        if isfield(params, 'allow_gamma_neg') && params.allow_gamma_neg
            params.gamma = clip(params.gamma, -1, 1);
        else
            params.gamma = clip(params.gamma, 0, 1);
        end
    elseif startsWith(fields{iF}, 'var_x', 'IgnoreCase', true)
        params.(fields{iF}) = clip(params.(fields{iF}), 1e-9, inf);
    elseif startsWith(fields{iF}, 'updates', 'IgnoreCase', true)
        params.(fields{iF}) = clip(round(params.(fields{iF})), 1, inf);
    elseif startsWith(fields{iF}, 'samples', 'IgnoreCase', true)
        params.(fields{iF}) = clip(round(params.(fields{iF})), 1, inf);
    elseif startsWith(fields{iF}, 'bound', 'IgnoreCase', true)
        params.(fields{iF}) = clip(params.(fields{iF}), 1e-9, inf);
    end
end

% During fitting of IS model, var_s and samples are highly correlated, so we reparameterize var_s as
% variance per sample.
if isfield(params, 'var_s_per_sample')
    params.var_s_per_sample = clip(params.var_s_per_sample, 1e-9, inf);
    params.var_s = params.var_s_per_sample * params.samples;
elseif isfield(params, 'var_s')
    params.var_s = clip(params.var_s, 1e-9, inf);
    params.var_s_per_sample = params.var_s / params.samples;
end
end

function val = clip(val, lo, hi)
val = min(max(val, lo), hi);
end