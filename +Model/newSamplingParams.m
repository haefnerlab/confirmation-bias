function params = newSamplingParams(varargin)
%NEWSAMPLINGPARAMS construct params struct for the sampling model.
%
% params = NEWSAMPLINGPARAMS('key1', val1, 'key2', val2, ...) create params
% with defaults and override any specified 'key' with the given value. For
% example, 'params = NEWSAMPLINGPARAMS('samples', 10)' creates a struct
% with all default params (e.g. params.var_x is 0.5), but with
% params.samples set to 10.
%
% Param keys are
%   var_e   - variance of gaussian p(e|x)
%   var_x   - variance of gaussian(s) p(x|D)
%   p_x     - weight of x modes, i.e. p(x|D) =
%             p_x*N(D,var_x)+(1-p_x)*N(-D,var_x)
%   prior_D - prior probability of D=+1
%   gamma   - value in [0,1]; how much the prior is 'subtracted
%             out'. 0 = low variance, 1 = unbiased.
%   samples - num samples per piece of evidence

params = struct(...
    'var_e', 0.1, ...
    'var_x', 0.5, ...
    'p_x', 1, ...
    'prior_D', 0.5, ...
    'gamma', 0, ...
    'samples', 1);

% Parse any extra (..., 'key', value, ...) pairs passed in through
% varargin.
for val_idx=2:2:length(varargin)
    key = varargin{val_idx-1};
    if ~ischar(key)
        warning('invalid input to newSamplingParams. All arguments should be (..., ''key'', value, ...)');
    elseif ~isfield(params, key)
        warning('unrecognized sampling parameter: ''%s''', key);
    else
        params.(key) = varargin{val_idx};
    end
end
end