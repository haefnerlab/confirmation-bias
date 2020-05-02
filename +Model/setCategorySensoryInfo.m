function params = setCategorySensoryInfo(params, ci, si)
% Set data-generating params
params.category_info = ci;
params.sensory_info = si;
% Set inference params
params.p_match = ci;
params.var_s = Model.getEvidenceVariance(si);
end