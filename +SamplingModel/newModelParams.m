function params = newModelParams(varargin)
%SAMPLINGMODEL.NEWMODELPARAMS construct params struct for the sampling model.
%
% params = SAMPLINGMODEL.NEWMODELPARAMS('key1', val1, 'key2', val2, ...) create params
% with defaults and override any specified 'key' with the given value. For
% example, 'params = SAMPLINGMODEL.NEWMODELPARAMS('samples', 10)' creates a struct
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
%   sensory_info  - probability of recovering x from e (see SamplingModel.getEvidenceVariance)
%   seed          - seed for generating data
% SAMPLING MODEL
%   var_s   - variance of gaussian p(s|x)
%   var_x   - variance of gaussian(s) p(x|C)
%   p_match - weight of x modes, i.e. p(x|C) = p_match*N(C,var_x)+(1-p_match)*N(-C,var_x)
%   prior_C - prior probability of C=+1
%   gamma   - value in [0,1]; how much the bias is 'subtracted out'.
%   samples - num "effective" samples of x per frame of evidence. Also the number of updates to p(C)
%             per frame.
%   batch   - effective "dimensionality" of x; how many independent measurements we get of x 
%             "per sample"
%  importance_norm - if true, uses normalized importance sampling for each batch. If false, uses
%                    'unbiased' weights.
%   noise   - amount of noise added, each frame, to accumulated log probabilities so far

params = struct(...
    'save_dir', fullfile('+SamplingModel', 'saved results'), ...
    'trials', 1000, ...
    'frames', 10, ...
    'category_info', .8, ...
    'sensory_info', .8, ...
    'seed', randi(1000000000), ...
    'var_s', 0.1, ...
    'var_x', 0.1, ...
    'p_match', 1, ...
    'prior_C', 0.5, ...
    'gamma', 0, ...
    'samples', 1, ...
    'batch', 1, ...
    'importance_norm', true, ...
    'noise', 0);

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