function params = setCategorySensoryInfo(params, ci, si)
% Set data-generating params
params.category_info = ci;
params.sensory_info = si;

% Set inference params
params.p_match = ci;
params.var_s = Model.getEvidenceVariance(si);

% In some special experiments we set 'gamma' as a function of the location in the space
if ~isempty(params.gammafun)
    params.gamma = params.gammafun(ci,si);
end
end