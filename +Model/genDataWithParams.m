function data = genDataWithParams(params, negtrials)
%Model.GENDATAWITHPARAMS generates a set of 'trials' (each is a 1xframes vector of real
%numbers), all with correct choice +1, with statistics matching the given sampling params.

if ~isempty(params.seed)
    rng(params.seed, 'twister');
end

% Generate the 'center' of each frame according to 'p_match'; with probability 'p_match' it is +1
% and with probability '1-p_match' it is -1.
centers = sign(params.category_info - rand(params.trials, params.frames));

% Use var_s as the variance of data around each center.
var_s = Model.getEvidenceVariance(params.sensory_info);

% Draw signal from around the center with stdev calculated above.
data = centers + randn(params.trials, params.frames) * sqrt(var_s);

if nargin >= 2 && negtrials
    flip = rand(params.trials, 1) < 0.5;
    data(flip, :) = -data(flip, :);
end
end