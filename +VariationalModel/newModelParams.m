function params = newModelParams(varargin)
%NEWMODELPARAMS construct params struct for the variational model.
%
% params = NEWMODELPARAMS('key1', val1, 'key2', val2, ...) create params with defaults and override
% any specified 'key' with the given value. For example, 'params = NEWMODELPARAMS('noise', .4)'
% creates a struct with all default params but with params.noise set to .4.
%
% Param keys are
% GENERAL CONFIGURATION
%   save_dir      - directory to save/load precomputed results
% DATA GENERATION
%   trials        - number of trials
%   frames        - number of frames per trial
%   category_info - probability that x matches c (or AROC of C from x in the gaussian model, determining var_x)
%   sensory_info  - probability of recovering x from s (see Model.getEvidenceVariance)
%   seed          - seed for generating data
% VARIATIONAL MODEL
%   model_fun - either @VariationalModel.runLatentZNormalX or ...Factorized
%   var_s   - variance of gaussian p(s|x)
%   var_x   - variance of gaussian(s) p(x|C)
%   prior_C - prior probability of C=+1
%   gamma   - value in [0,1]; how much the bias is 'subtracted out'.
%   updates - number of 'updates' per frame
%   noise   - amount of noise added, each frame, to accumulated log probabilities so far

params = struct(...
    'save_dir', fullfile('+VariationalModel', 'saved results'), ...
    'trials', 1000, ...
    'frames', 10, ...
    'category_info', .8, ...
    'sensory_info', .8, ...
    'seed', randi(1000000000), ...
    'model_fun', @VariationalModel.runLatentZNormalX, ...
    'var_s', .1, ...
    'var_x', .1, ...
    'p_match', .8, ...
    'prior_C', .5, ...
    'gamma', 0, ...
    'updates', 10, ...
    'noise', 0);

% Parse any extra (..., 'key', value, ...) pairs passed in through varargin.
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