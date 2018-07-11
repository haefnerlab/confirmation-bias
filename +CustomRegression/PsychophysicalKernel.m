function [weights, postVal, errors, map_ridge, map_ar1, map_curvature] = ...
    PsychophysicalKernel(data, responses, hpr_ridge, hpr_ar1, hpr_curvature, standardize)
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
% ... = PSYCHOPHYSICALKERNEL(..., standardize) sets standardization rule.
% When standardize == 0, no preprocessing is done. When standardize == 1,
% the data are standardized assuming iid and zero-mean. When standardize ==
% 2, the data are standardized using the built-in zscore() function.
% Defaults to 0.

if nargin < 6, standardize = 0; end
if nargin < 5, hpr_curvature = 0; end
if nargin < 4, hpr_ar1 = 0; end
if nargin < 3, hpr_ridge = 0; end

% Standardize each regressor.
switch standardize
    case 0
        % do nothing
    case 1
        % assume 0 mean (nothing to subtact) and iid (std taken over all data)
        data = data / std(data(:));
    case 2
        data = zscore(data);
    otherwise
        error('Expected argument ''standardize'' to be one of [0, 1, 2]');
end

% Add a column of ones to data for a bias term.
data = [data ones(size(data,1),1)];

% convert boolean to float type
assert(islogical(responses));
responses = 1.0 * responses(:);

% Grid search will be done over all combinations (i.e. the cartesian
% product) of given hyperparameters.
grid_size = [length(hpr_ridge) length(hpr_ar1) length(hpr_curvature)];
n_gridpts = prod(grid_size);

% Each entry in 'results' will itself be a cell array containing the return
% values (i.e. {weights, postval, errors, ...})
results = cell(n_gridpts, 1);

for i=1:n_gridpts
    % Determine which hyperparameters to use this iteration by treating i
    % as a 1d index into the 3d grid of ridge/ar1/curvature values.
    [idx_ridge, idx_ar1, idx_curvature] = ind2sub(grid_size, i);
    
    [weights, postVal, errors] = do_fit(data, responses, hpr_ridge(idx_ridge), hpr_ar1(idx_ar1), hpr_curvature(idx_curvature));
    
    % Record all results for this set of hyperparameters.
    results{i} = {weights, postVal, errors, hpr_ridge(idx_ridge), hpr_ar1(idx_ar1), hpr_curvature(idx_curvature)};
end

% Find and return the MAP result.
postVals = cellfun(@(result) result{2}, results);
[~, idx_map] = max(postVals);
[weights, postVal, errors, map_ridge, map_ar1, map_curvature] = results{idx_map}{:};

end

function [optim_weights, postVal, errors] = do_fit(data, responses, hpr_ridge, hpr_ar1, hpr_curvature)
% Perform fitting at a single setting of hyperparameters / at a single
% point in the grid search.

[trials, p] = size(data);

break_points = p-1; % don't smooth onto bias term
D1 = derivative_matrix(p, break_points);
D2 = second_derivative_matrix(p, break_points);
prior_matrix = hpr_ridge * eye(p) + hpr_ar1 * (D1'*D1) + hpr_curvature * (D2'*D2);

    function [nlp, grad, hessian] = neg_log_posterior(weights)
        % COPIED FROM CustomRegression.PsychophysicalKernel
        neg_log_prior = 0.5 * weights' * prior_matrix * weights;
        nlp = neg_log_prior + neg_bernoulli_log_likelihood(data, responses, weights);
        
        if nargout >= 2
            logits = data * weights;
            ps = exp(logits) ./ (1 + exp(logits));
            grad = prior_matrix * weights + data' * ps - data' * responses;
        end
        
        if nargout >= 3
            % get sigmoid derivative at each trial
            deriv_ps = ps .* (1 - ps);
            % Using sparse diagonal matrix to avoid creating trials x
            % trials full diagonal matrix.
            sp_diag_deriv_ps = spdiags(deriv_ps, 0, spalloc(trials, trials, trials));
            hessian = prior_matrix + data' * sp_diag_deriv_ps * data;
        end
    end

% Fit weights using 'fminunc'
options = optimoptions('fminunc', ...
    'Display', 'none', ...
    'Algorithm', 'trust-region', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');
[optim_weights, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_log_posterior, zeros(p, 1), options);
postVal = -negPostVal;
% attempt to invert the hessian for standard error estimate - this sometimes fails silently,
% returning NaN.
errors = sqrt(diag(inv(hessian)));
end

function D = derivative_matrix(n, break_at)
D = full(CustomRegression.create_forward_difference_matrix([1 n], [-1 1]));
for i=break_at
    D(i,:) = 0;
end
end

function D = second_derivative_matrix(n, break_at)
D = full(CustomRegression.create_forward_difference_matrix([1 n], [1 -2 1]));
for i=break_at
    D(i-1, :) = 0;
    D(i, :) = 0;
end
end

function LL = neg_bernoulli_log_likelihood(data, responses, weights)
logits = data * weights;
% we want \sum_i responses(i)*log(p(i)) + (1-responses(i))*log(1-p(i))
% where p(i) = (1+exp(-logit(i)))^-1. However, this can lead to numerical
% instability due to log(exp(...)). Reparameterizing, note that p^r is 1
% when r=0 and exp(logit)/(1+exp(logit)) when r=1, and likewise (1-p)^(1-r)
% is simply 1 when r=1 and 1/(1+exp(logit)) when r=0. Hence p^r (1-p)^(1-r)
% becomes exp(r*logit)/(1+exp(logit)), and negative log of this is simply:
log_bernoulli = -responses .* logits(:) + log(1 + exp(logits(:)));
LL = sum(log_bernoulli);
end
