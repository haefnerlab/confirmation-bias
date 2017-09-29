function params = newModelParams(varargin)
%NEWMODELPARAMS construct params struct for the sampling model.
%
% params = NEWMODELPARAMS('key1', val1, 'key2', val2, ...) create params
% with defaults and override any specified 'key' with the given value. For
% example, 'params = NEWMODELPARAMS('samples', 10)' creates a struct
% with all default params (e.g. params.var_x is 0.5), but with
% params.samples set to 10.
%
% Param keys are
% GENERAL CONFIGURATION
%   save_dir      - directory to save/load precomputed results
% DATA GENERATION
%   trials        - number of trials
%   frames        - number of frames per trial
%   category_info - probability that x matches c
%   sensory_info  - probability of recovering x from e (see Model.getEvidenceVariance)
%   seed          - seed for generating data
% SAMPLING MODEL
%   var_e   - variance of gaussian p(e|x)
%   var_x   - variance of gaussian(s) p(x|C)
%   p_match - weight of x modes, i.e. p(x|C) =
%             p_match*N(C,var_x)+(1-p_match)*N(-C,var_x)
%   prior_C - prior probability of C=+1
%   gamma   - value in [0,1]; how much the prior is 'subtracted
%             out'. 0 = low variance, 1 = unbiased.
%   samples - num "effective" samples of x per frame of evidence. Also the
%             number of updates to p(C) per frame.
%   batch   - effective "dimensionality" of x; how many independent 
%             measurements we get of x "per sample" 

params = struct(...
    'save_dir', fullfile('+Model', 'saved results'), ...
    'trials', 1000, ...
    'frames', 10, ...
    'category_info', .8, ...
    'sensory_info', .8, ...
    'seed', randi(1000000000), ...
    'var_e', 0.1, ...
    'var_x', 0.5, ...
    'p_match', 1, ...
    'prior_C', 0.5, ...
    'gamma', 0, ...
    'samples', 1, ...
    'batch', 1);

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