function [data, descriptor] = genDataWithParams(trials, frames, params)
%GENDATAWITHPRIORLIKELIHOOD generates a set of 'trials' (a 1xframes vector
%of real numbers), all with correct choice +1, with statistics matching the
%given sampling params.

% generate the 'center' of each frame according to 'p_match'; with
% probability 'p_match' it is +1 and with probability '1-p_match' it is -1.
centers = sign(params.p_match - rand(trials, frames));

% The net variance of p(e|D) is var_e + var_x.
stdev = sqrt(params.var_x + params.var_e);

% draw signal from around the center with stdev calculated above.
data = centers + randn(trials, frames) * stdev;

descriptor = sprintf('%dx%d_all_params', trials, frames);
end