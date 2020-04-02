function tf = isStochastic(params, noisefloor)
%MODEL.ISSTOCHASTIC return whether the given model params are a stochastic integrator, i.e. will
%output different log odds on multiple runs given the same stimuli.
if nargin < 2, noisefloor = 1e-9; end
tf = strcmpi(params(1).model, 'is') || ...
    (any(strcmpi(params(1).model, {'vb', 'vb-czx', 'itb'})) && params(1).noise > noisefloor);
end