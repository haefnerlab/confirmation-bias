function [weights, postVal, errors, map_ridge, map_ar1, map_curvature] = ...
    PsychophysicalKernel(data, responses, hpr_ridge, hpr_ar1, hpr_curvature, split_smoothness)
%PSYCHOPHYSICALKERNEL Regress PK
%
% [ weights, postVal, errors ] = PSYCHOPHYSICALKERNEL(data, responses) 
% where data is [trials x regressors] and responses is [trials x 1] of
% booleans, computes regression weights 'weights', the value of the
% log posterior, and the estimated errors on each weight. Errors may be NaN
% if the Hessian is degenerate at the MAP solution.
%
% PSYCHOPHYSICALKERNEL(data, responses, hpr_ridge, hpr_ar1, hpr_curvature)
% performs grid search over all combinations of the given prior
% hyperparameters (arrays of values). Ridge controls the magnitude of the
% weights (a bias towards zero, i.e. ridgre regression). AR1 is a prior on
% the first temporal derivative of the weights, encouraging them to be
% close to constant. Curvature is a prior on the second temporal derivative
% of the weights, encouraging them to be smooth. For example,
% PSYCHOPHYSICALKERNEL(D, R, [0 0.5 1], [0], [10 100]) will fit using all
% of [0 0 10], [0 0 100], [.5 0 10], [.5 0 100], [1 0 10], [1 0 100] for
% the ridge/ar1/curvature hyperparameters respectively. Three additional
% return values contain the MAP estimate of each hyperparameter..
% [ w, p, e, map_ridge, map_ar1, map_curvature ] = PSYCHOPHYSICALKERNEL(...)
%
% PSYCHOPHYSICALKERNEL(..., split_smoothness) when split_smoothness is
% true, the smoothness constraints are not enforced at the midpoint of the
% data. This is useful when computing temporal weights for two indpendent
% 'sensors' with concatenated data, e.g. audio data to the left and right
% ear separately.

if nargin < 6, split_smoothness = false; end
if nargin < 5, hpr_curvature = 0; end
if nargin < 4, hpr_ar1 = 0; end
if nargin < 3, hpr_ridge = 0; end

% convert boolean to float type
responses = 1.0 * responses;

[~, p] = size(data);

break_points = p-1; % don't smooth onto bias term
if split_smoothness, break_points = [break_points floor(p/2)]; end
D1 = derivative_matrix(p, break_points);
D2 = second_derivative_matrix(p, break_points);

% Grid search will be done over all combinations (i.e. the cartesian
% product) of given hyperparameters.
grid_size = [length(hpr_ridge) length(hpr_ar1) length(hpr_curvature)];
n_gridpts = prod(grid_size);

% Each entry in 'results' will itself be a cell array containing the return
% values (i.e. {weights, postval, errors, ...})
results = cell(n_gridpts, 1);

compute_error = nargout > 2;

parfor i=1:n_gridpts
    % Determine which hyperparameters to use this iteration by treating i
    % as a 1d index into the 3d grid of ridge/ar1/curvature values.
    [idx_ridge, idx_ar1, idx_curvature] = ind2sub(grid_size, i);
    % Construct a loss function given these hyperparameters.
    negLogPost = @(w) ...
        - log_prior(w, D1, D2, hpr_ridge(idx_ridge), hpr_ar1(idx_ar1), hpr_curvature(idx_curvature)) ...
        - bernoulli_log_likelihood(data, responses, w);

    % Fit weights using 'fminunc', only computing the hessian (which is
    % slow) if errors are requested.
    if compute_error
        [weights, negPostVal, ~, ~, ~, hessian] = fminunc(negLogPost, zeros(p, 1));
        % attempt to invert the hessian for error bars - this sometimes
        % fails silently, returning NaN.
        errors = abs(diag(inv(hessian)));
    else
        [weights, negPostVal] = fminunc(negLogPost, zeros(p, 1));
        errors = zeros(size(weights));
    end
    postVal = -negPostVal;
    
    % Record all results for this set of hyperparameters.
    results{i} = {weights, postVal, errors, hpr_ridge(idx_ridge), hpr_ar1(idx_ar1), hpr_curvature(idx_curvature)}
end

% Find and return the MAP result.
postVals = cellfun(@(result) result{2}, results);
[~, idx_map] = max(postVals);
[weights, postVal, errors, map_ridge, map_ar1, map_curvature] = results{idx_map}{:};

end

function D = derivative_matrix(n, break_at)
D = - eye(n) + diag(ones(n-1, 1), 1);
D(end) = 0;
for i=break_at
    D(i,:) = 0;
end
end

function D = second_derivative_matrix(n, break_at)
D = eye(n) - 2*diag(ones(n-1, 1), 1) + diag(ones(n-2, 1), 2);
D(end) = 0;
D(end-1, end-2:end) = 0;
for i=break_at
    D(i-1, :) = 0;
    D(i, :) = 0;
end
end

function LP = log_prior(weights, D1, D2, hpr_ridge, hpr_ar1, hpr_curvature)
ridge = -0.5 * dot(weights, weights);
ar1 = 0;
if hpr_ar1 > 0
    ar1 = -0.5 * dot(D1*weights, D1*weights);
end
curvature = 0;
if hpr_curvature > 0
    curvature = -0.5 * dot(D2*weights, D2*weights);
end
LP = hpr_ridge * ridge + hpr_ar1 * ar1 + hpr_curvature * curvature;
end

function LL = bernoulli_log_likelihood(data, responses, weights)
p = sigmoid(data * weights);
LL = dot(responses, log(p)) + dot(1-responses, log(1-p));
end

function y = sigmoid(x)
y = (1 + exp(-x)).^-1;
end