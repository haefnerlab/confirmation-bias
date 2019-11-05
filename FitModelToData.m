% function [bestfit_params, loglike] = FitModelToData(SubjectData, varargin)

init_params = Model.newModelParams('model', 'is', ...
    'gamma', 0.1, ...
    'samples', 5, ...
    'updates', 2, ...
    'lapse', 0.01);

nInner = 10;

sensor_noise_grid = logspace(-2, 1, 11);

for iNoise=length(sensor_noise_grid):-1:1
    [param_set, stim_set, choice_set] = SubjectDataToModelParams(SubjectData, ...
        sensor_noise_grid(iNoise), init_params);
    
    fields = {'prior_C', 'gamma', 'samples', 'lapse'};
    LB = [0 0 1 0];
    UB = [1 1 100 1];
    PLB = LB;
    PUB = UB;
    
    x0 = cellfun(@(f) init_params.(f), fields);
    
    wrapper_args = {@Fitting.choiceModelLogLikelihood, param_set, fields, {stim_set, choice_set, nInner}};
    
    vbmc_options = vbmc('defaults');
    vbmc_options.UncertaintyHandling = 'yes';
    [VP(iNoise), ELBO(iNoise), ELBO_SD(iNoise), EXITFLAG(iNoise)] = vbmc(@loglikefn_wrapper, ...
        x0, LB, UB, PLB, PUB, vbmc_options, wrapper_args{:});
    
end

% end

function [val] = loglikefn_wrapper(xval, loglikefn, params, fields, args)
for j=1:length(params)
    params(j).var_s_per_sample = params(j).var_s / params(j).samples;
    for i=1:length(fields)
        params(j).(fields{i}) = xval(i);
    end
end
val = loglikefn(params, args{:});
end