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
        params.gamma = clip(params.gamma, -1, 1);
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
end

function val = clip(val, lo, hi)
val = min(max(val, lo), hi);
end