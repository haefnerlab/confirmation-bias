function params = newModelParams(varargin)
%MODEL.NEWMODELPARAMS construct params struct for the sampling model.
%
% params = MODEL.NEWMODELPARAMS('key1', val1, 'key2', val2, ...) create params with defaults and
% override any specified 'key' with the given value. For example, 
%'params = MODEL.NEWMODELPARAMS('samples', 10)' creates a struct with all default params (e.g. 
% params.gamma is 0), but with params.samples set to 10.
%
% Param keys are
% GENERAL CONFIGURATION
%   save_dir      - directory to save/load precomputed results
%   model         - one of 'is', 'vb', 'vb-czx', or 'ideal' specifying which model to use
% DATA GENERATION
%   trials        - number of trials
%   frames        - number of frames per trial
%   category_info - probability that x matches c
%   sensory_info  - probability of recovering x from e (see Model.getEvidenceVariance)
%   seed          - seed for generating data
% MODEL PARAMETERS
%   var_s   - variance of gaussian p(s|x)
%   p_match - weight of x modes, i.e. p(x|C) = p_match*N(C,var_x)+(1-p_match)*N(-C,var_x)
%   prior_C - prior probability of C=+1
%   var_x   - variance of each mode of p(x|C)
%   gamma   - value in [0,1]; how much the bias is 'subtracted out'.
%   updates - the number of updates to p(C) per frame, corresponding to the brain's sampling time.
%   noise   - amount of noise added, each frame, to accumulated log probabilities so far
%   lapse   - lapse rate (percent of trials making a random choice regardless of stimulus)
%   samples - effective "dimensionality" of x; how many independent measurements we get of x 
%             per update (e.g. multiple independent sampling chains). Only used in IS model.
%   step_size - What percent of the "full update" to log(q(C=+1)/q(C=-1)) to apply each step. Only
%               used in VB model.

params = struct(...
    'save_dir', fullfile('+Model', 'saved results'), ...
    'model', 'is', ...
    'trials', 10000, ...
    'frames', 10, ...
    'category_info', .75, ...
    'sensory_info', .75, ...
    'seed', randi(1000000000), ...
    'var_s', Model.getEvidenceVariance(.75), ...
    'p_match', .75, ...
    'prior_C', 0.5, ...
    'var_x', 0.1, ...
    'gamma', 0, ...
    'updates', 5, ...
    'noise', 0, ...
    'lapse', 0, ...
    'samples', 5, ...
    'step_size', 0.1);

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

%% Sanity check

if strcmp(params.model, 'ideal') && params.updates > 1
    warning('Ideal observer; setting params.updates to 1');
    params.updates = 1;
end
end