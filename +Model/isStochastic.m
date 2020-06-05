function tf = isStochastic(params, noisefloor)
%MODEL.ISSTOCHASTIC return whether the given model params are a stochastic integrator, i.e. will
%output different log odds on multiple runs given the same stimuli.
if nargin < 2, noisefloor = 1e-9; end
assert(ismember(params(1).model, {'is', 'vb', 'vb-czx', 'itb', 'itb-int', 'ideal'}));
% IS model is inherently stochastic. Other models are only stochastic if noise > 0. The itb-int
% model is deterministic in the sense that it always (approximately) integrates out the noise in
% exactly the same way.
tf = strcmpi(params(1).model, 'is') || ...
    (any(strcmpi(params(1).model, {'vb', 'vb-czx', 'itb'})) && params(1).noise > noisefloor);
end