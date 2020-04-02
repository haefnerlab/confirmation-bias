function name = getModelStringID(params, drop_cs_terms)
switch lower(params.model)
    case 'is'
        name = sprintf('is_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_ns%d_noise%.2e_temp%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, params.samples, ...
            params.noise, params.temperature, params.lapse);
        if ~params.importance_norm
            name = [name '_unnorm'];
        end
    case 'vb'
        name = sprintf('vb_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_step%.2e_noise%.2e_temp%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, ...
            params.step_size, params.noise, params.temperature, params.lapse);
    case 'vb-czx'
        name = sprintf('vb_cxz_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_step%.2e_noise%.2e_temp%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, ...
            params.step_size, params.noise, params.temperature, params.lapse);
    case 'ideal'
        name = sprintf('ideal_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_temp%.2e_lapse%.2e', params.trials, ...
            params.frames, params.category_info, params.sensory_info, params.var_s, params.p_match, ...
            params.prior_C, params.var_x,  params.temperature, params.lapse);
    case 'itb'
        name = sprintf('itb_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_bnd%.2f_gam%.2f_noise%.2e_temp%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.bound, params.gamma, params.noise, ...
            params.temperature, params.lapse);
    otherwise
        error('Unrecognized model type: %s', params.model);
end

if nargin >= 2 && drop_cs_terms
    % Drop terms corresponding to individual points in the space..
    name = regexprep(name, '_cinfo[0-9.]+_sinfo[0-9.]+_vs[0-9.]+_pm[0-9.]+', '');
    
    % But maybe add a term describing how gamma changes over the space
    if isfield(params, 'gamma_min') && ~isempty(params.gamma_min)
        name = regexprep(name, '_gam[0-9.]+', sprintf('_gmin%.2f_gmax%.2f', params.gamma_min ,params.gamma_max));
    end
end

end