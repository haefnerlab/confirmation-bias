function name = getModelStringID(params, ideal_observer)
if nargin < 2, ideal_observer = false; end
if ~ideal_observer
    if params.importance_norm
        norm_str = '_norm';
    else
        norm_str = '_unnorm';
    end
    name = sprintf('%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_ns%d_nb%d_noise%.2e%s', ...
        params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
        params.p_match, params.prior_C, params.var_x, params.gamma, params.samples, params.batch, ...
        params.noise, norm_str);
else
    name = sprintf('%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f', ...
        params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
        params.p_match, params.prior_C);
end
end