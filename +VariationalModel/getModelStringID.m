function name = getModelStringID(params)
model_fun_info = functions(params.model_fun);
fun_name_parts = split(model_fun_info.function, '.');
fun_name = fun_name_parts{end};
name = sprintf('%dx%d_cinfo%.3f_sinfo%.3f_%s_vs%.2f_vx%.2f_pm%.2f_pC%.2f_u%d_gam%.2f_nz%.2x', ...
    params.trials, params.frames, params.category_info, params.sensory_info, fun_name, params.var_s, ...
    params.var_x, params.p_match, params.prior_C, params.updates, params.gamma, params.noise);
end