function [gt_params, sigs, choices] = GetGroundTruthSimData(model, phases)
switch upper(model)
    % Note: seeds make result repeatable even if not cached, since calling 'runVectorized'
    % immediately after 'genDataWithParams' below causes 'runVectorized' to behave according to the
    % same rng seed.
    case 'IS'
        base_params = Model.newModelParams('model', 'is', 'trials', 1000, 'temperature', 0.1, 'var_x', 0.1, 'gamma', 0.1, 'samples', 5, 'updates', 5, 'seed', 18863457);
    case 'ITB'
        base_params = Model.newModelParams('model', 'itb', 'trials', 1000, 'temperature', 0.1, 'bound', 1.2, 'gammafun', @(ci,si) 1-ci, 'noise', 0.35, 'updates', 1, 'seed', 18863457);
    otherwise
        error('Bad name');
end

assert(all(ismember(phases, [1 2])), 'Bad phase identifier Use 1 for HSLC and 2 for LSHC!');

sens_cat_pts = Model.getThresholdPoints(0.51:0.02:0.99, base_params, .7, 5);

% HSLC ground truth
gt_params(1) = Model.setCategorySensoryInfo(base_params, sens_cat_pts(end,2), sens_cat_pts(end,1));
gt_data{1} = Model.genDataWithParams(gt_params(1));
gt_sim(1) = Model.runVectorized(gt_params(1), gt_data{1});

% LSHC ground truth
gt_params(2) = Model.setCategorySensoryInfo(base_params, sens_cat_pts(1,2), sens_cat_pts(1,1));
gt_data{2} = Model.genDataWithParams(gt_params(2));
gt_sim(2) = Model.runVectorized(gt_params(2), gt_data{2});

% Select by phase.. kind of hackily using the flag for HSLC vs LSHC (1 vs 2 respectively) as
% indices.
if length(phases) == 1
    % Select a signle phase
    gt_params = gt_params(phases);      % struct
    sigs = gt_data{phases};             % plain old matrix
    choices = gt_sim(phases).choices;   % column vector
else
    % Select concatenated phases in some order (e.g. phases may be [1 2] or [2 1])
    gt_params = gt_params(phases);      % struct array
    sigs = gt_data(phases);             % cell array
    choices = {gt_sim(phases).choices}; % cell array
end
end