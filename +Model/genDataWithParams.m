function [data, categories,seed] = genDataWithParams(params)
%Model.GENDATAWITHPARAMS generates a set of 'trials' (each is a 1xframes vector of real
%numbers), with statistics matching the given sampling params.

if isfield(params, 'seed') && ~isempty(params.seed)
	seed = params.seed;
else
	seed = randi(1e9);
end
rng(seed, 'twister');

% Generate the 'center' of each frame according to 'p_match'; with probability 'p_match' it is +1
% and with probability '1-p_match' it is -1.
centers = sign(params.category_info - rand(params.trials, params.frames));

% Use var_s as the variance of data around each center.
var_s = Model.getEvidenceVariance(params.sensory_info);

% Draw signal from around the center with stdev calculated above.
data = centers + randn(params.trials, params.frames) * sqrt(var_s);

% So far, all trials have 'correct' category +1. Randomly flip half of them.
categories = ones(params.trials, 1);
flip = rand(params.trials, 1) < 0.5;
data(flip, :) = -data(flip, :);
categories(flip) = -1;
end