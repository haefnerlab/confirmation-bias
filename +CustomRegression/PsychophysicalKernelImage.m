function [ weights, postVal, errors ] = PsychophysicalKernelImage(data, responses, ...
    hpr_ridge, hpr_ar1, hpr_curvature, ...
    hpr_sp_ridge, hpr_sp_ar1, hpr_sp_curvature, ...
    init_template)
%PSYCHOPHYSICALKERNELIMAGE Joint spatial and temporal PsychophysicalKernel.
%
% [weights, postVal, errors] = PSYCHOPHYSICALKERNELIMAGE(data, responses)
% Computes separate spatial/temporal weights, the posterior value, and an
% estimate of errors on each weight. 'data' is a cell array of images per
% 'frame' (with shape [height width frames]) and responses is a boolean
% vector of the subject's choices. Each image must be the same size. If a
% trial contains images with P pixels over T frames, then 'weights' will
% have length P+T such that 'reshape(weights(1:P), [height width])' is the
% spatial kernel, weights(P+1:end-1) is the temporal kernel, and
% weights(end) is a bias term.
%
% Temporal hyperparameters and grid search are the same as in
% PsychophysicalKernel(). All unspecified hyperparameters default to 0.
% [weights, postVal, errors, map_ridge, map_ar1, map_curvature] = ...
%   PSYCHOPHYSICALKERNELIMAGE(data, responses, hpr_ridge, hpr_ar1, hpr_curvature)
%
% Additional analogous spatial hyperparameters are available as well:
% [~, ~, ~, ~, ~, ~, map_sp_ridge, map_sp_ar1, map_sp_curvature] = ...
%   PSYCHOPHYSICALKERNELIMAGE(data, responses, hpr_ridge, hpr_ar1, hpr_curvature, ...
%   hpr_sp_ridge, hpr_sp_ar1, hpr_sp_curvature)
%
% these hpr_sp_* hyperparameters control put a small-magnitude prior on the
% spatial weights, and their first and second spatial derivatives (in both
% x and y directions) respectively.
%
% PSYCHOPHYSICALKERNELIMAGE(..., init_template) initializes the spatial
% template to the given image (which must have the same size as the images
% in the data)
[height, width, frames] = size(data{1});
trials = length(data);

% Standardize each regressor.
data = cellfun(@(trial) reshape(trial, [1 height*width*frames]), data, ...
    'UniformOutput', false);
data = zscore(cat(1, data{:}));
data = reshape(data, [trials height width frames]);
data = arrayfun(@(i) squeeze(data(i, :, :, :)), 1:trials, 'UniformOutput', false);

if nargin < 9, init_template = 0.01 * randn(height, width); end
if nargin < 8, hpr_sp_curvature = 0; end
if nargin < 7, hpr_sp_ar1 = 0; end
if nargin < 6, hpr_sp_ridge = 0; end
if nargin < 5, hpr_curvature = 0; end
if nargin < 4, hpr_ar1 = 0; end
if nargin < 3, hpr_ridge = 0; end

% convert boolean to float type
responses = 1.0 * responses;

pixels = height*width;
% Reshape each image into a column vector so that each trial is a (pixels x
% frames) matrix.
flat_data = cellfun(@(trial) reshape(trial, [pixels, frames]), data, 'UniformOutput', false);
init_weights = [init_template(:); 0.01 * randn(frames + 1, 1)];

% Construct temporal priors.
break_points = frames-1; % don't smooth onto bias term
Dt = derivative_matrix(frames, break_points);
Dtt = second_derivative_matrix(frames, break_points);

% Construct spatial priors.
[Dx, Dy] = spatial_derivative_matrix(height, width);
[Dxx, Dyy, Dxy] = spatial_second_derivative_matrix(height, width);

% Create loss function (negative log posterior)
function temporal_signal = dot_spatial(weights)
    projected_cell = cellfun(@(trial) weights'*trial, flat_data, 'UniformOutput', false);
    temporal_signal = [vertcat(projected_cell{:}) ones(length(flat_data), 1)];
end

negLogPost = @(w) ...
    - log_prior(w(pixels+1:end-1), Dt, Dtt, hpr_ridge, hpr_ar1, hpr_curvature) ... % temporal smoothness priors
    - log_spatial_prior(w(1:pixels), Dx, Dy, Dxx, Dyy, Dxy, hpr_sp_ridge, hpr_sp_ar1, hpr_sp_curvature) ... spatial smoothness priors
    - bernoulli_log_likelihood(dot_spatial(w(1:pixels)), responses, w(pixels+1:end)); % bernoulli likelihood, passing data first through spatial then through temporal kernels.

if nargout > 2
    [weights, negPostVal, ~, ~, ~, hessian] = fminunc(negLogPost,init_weights);
    % attempt to invert the hessian for error bars
    errors = abs(diag(inv(hessian)));
else
    [weights, negPostVal] = fminunc(negLogPost, init_weights);
end

postVal = -negPostVal;

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

function [Dx, Dy] = spatial_derivative_matrix(height, width)
Dx = CustomRegression.create_forward_difference_matrix([height width], [-1 1]);
Dy = CustomRegression.create_forward_difference_matrix([height width], [-1 1]');
end

function [Dxx, Dyy, Dxy] = spatial_second_derivative_matrix(height, width)
Dxx = CustomRegression.create_forward_difference_matrix([height width], [1 -2 1]);
Dyy = CustomRegression.create_forward_difference_matrix([height width], [1 -2 1]');
Dxy = CustomRegression.create_forward_difference_matrix([height width], [1 -1; -1 1]);
end

function LP = log_prior(weights, Dt, Dtt, hpr_ridge, hpr_ar1, hpr_curvature)
ridge = -0.5 * (sqrt(dot(weights, weights)) - 1)^2;
ar1 = 0;
if hpr_ar1 > 0
    ar1 = -0.5 * dot(Dt*weights, Dt*weights);
end
curvature = 0;
if hpr_curvature > 0
    curvature = -0.5 * dot(Dtt*weights, Dtt*weights);
end
LP = hpr_ridge * ridge + hpr_ar1 * ar1 + hpr_curvature * curvature;
end

function LP = log_spatial_prior(weights, Dx, Dy, Dxx, Dyy, Dxy, hpr_ridge, hpr_ar1, hpr_curvature)
ridge = -0.5 * dot(weights, weights);
ar1 = 0;
if hpr_ar1 > 0
    ar1 = -0.5 * (dot(Dx*weights, Dx*weights) + dot(Dy*weights, Dy*weights));
end
curvature = 0;
if hpr_curvature > 0
    curvature = -0.5 * (dot(Dxx*weights, Dxx*weights) + dot(Dyy*weights, Dyy*weights) + 2*dot(Dxy*weights, Dxy*weights));
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

