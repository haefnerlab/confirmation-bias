function [ varargout ] = PsychophysicalKernelGridSearch( data, responses, pkfun, varargin )
%PSYCHOPHYSICALKERNELGRIDSEARCH A wrapper around PsychophysicalKernel
%functions that does a grid search over hyperparameters in parallel.
%
% [W, LP, E, ...] = PSYCHOPHYSICALKERNELGRIDSEARCH(data, responses, pkfun, hprs1, hprs2, ...)
% fits the given pkfun (e.g. CustomRegression.PsychophysicalKernel) to the
% given data and responses, searching over the given values of
% hyperparameters. The number of grid points searched is
% length(hprs1)*length(hprs2)*... Output will be the outputs of pkfun, plus
% values of the hyperparameters at the MAP grid point. For example,
%
% fn = @CustomRegression.PsychophysicalKernel
% [weights, postval, errors, ridge, ar1, curvature] = ...
% PSYCHOPHYSICALKERNELGRIDSEARCH(data, responses, fn, ridges, ar1s, curvatures)
%
% Other non-grid-searched arguments to the pkfun can be passed in after a
% 'args' delimeter, e.g.
% 
% fn = @CustomRegression.PsychophysicalKernel
% PsychophysicalKernelGridSearch(d, r, fn, [0 .01 .1 1], [0], [0], 'args', true)
%
% ..will search over values [0, 0.01, 0.1 1] for hpr_ridge, and will set
% split_smoothness to true.

split_args = strcmpi(varargin, 'args');
extra_args = {};
if any(split_args)
    extra_args = varargin(find(split_args)+1:end);
    varargin = varargin(1:find(split_args)-1);
end

n_hprs = length(varargin);
all_hprs = cell(1, n_hprs);
[all_hprs{:}] = meshgrid(varargin{:});
all_hprs = cellfun(@(hpr_set) hpr_set(:), all_hprs, 'UniformOutput', false);
all_hprs = horzcat(all_hprs{:});
% all_hprs is now an N x H matrix, where N=num grid points and H=num
% hyperparameters. each row of all_hprs is a unique combination of the
% hyperparameters, i.e. a grid point.

n_gridpts = size(all_hprs,1);
fprintf('Beginning parallel search over %d gridpoints\n', n_gridpts);
results = cell(n_gridpts, 1);
parfor i=1:n_gridpts
    res = cell(1, 3+n_hprs);
    args = horzcat(num2cell(all_hprs(i,:)), extra_args);
    [weights, postVal, errors] = pkfun(data, responses, args{:});
    res{1} = weights;
    res{2} = postVal;
    res{3} = errors;
    [res{4:4+n_hprs-1}] = args{1:n_hprs};
    results{i} = res;
end

% Find MAP result
postVals = cellfun(@(result) result{2}, results);
[~, maxIdx] = max(postVals);
varargout = results{maxIdx};
end

