function [best_hprs, log_likelihoods] = xValidatePK(data, responses, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds)
%CUSTOMREGRESSION.XVALIDATEPK Perform folded cross-validation on logistic regression with varying
%hyperparameters. Evaluates the cartesian product of hpr_ridge x hpr_ar1 x hpr_curvature. Each term
%is evaluated over 'folds' splits.

% Begin by randomly shuffling both 'data' and 'responses' so that the x-validation splits are not
% sensitive to any ordering of the data.
[trials, ~] = size(data);
shuffle = randperm(trials);
data = data(shuffle, :);
responses = responses(shuffle);

% Standardize each regressor first, across the whole dataset (this is identical to the
% CustomRegression.PsychophysicalKernel function)
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

% Get start, end indices of each split
splitStart = round(linspace(1, trials+1, folds+1));
splitEnd = splitStart(2:end) - 1;

log_likelihoods = zeros(length(hpr_ridge), length(hpr_ar1), length(hpr_curvature), folds);
sz = size(log_likelihoods);

parfor ii=1:numel(log_likelihoods)
    [iRidge, iAR1, iCurve, iFold] = ind2sub(sz, ii);
    
    dataIdx = [1:splitStart(iFold)-1 splitEnd(iFold)+1:trials];
    weights = CustomRegression.PsychophysicalKernel(data(dataIdx, :), responses(dataIdx), ...
        hpr_ridge(iRidge), hpr_ar1(iAR1), hpr_curvature(iCurve), 0);
    
    log_likelihoods(ii) = bernoulli_log_likelihood(data(splitStart(iFold):splitEnd(iFold), :), ...
        responses(splitStart(iFold):splitEnd(iFold)), weights);
end

avg_ll = mean(log_likelihoods, 4);
[~, imax] = max(avg_ll(:));
[iRidge, iAR1, iCurve] = ind2sub(sz(1:3), imax);
% Err on the side of less regularization by choosing smoothing that is one order of magnitude less
% than the best.
iRidge = max(iRidge-1, 1);
iAR1 = max(iAR1-1, 1);
iCurve = max(iCurve-1, 1);
best_hprs = [hpr_ridge(iRidge), hpr_ar1(iAR1), hpr_curvature(iCurve)];

end

function LL = bernoulli_log_likelihood(data, responses, weights)
% Copied from of neg_bernoulli_log_likelihood() found in CustomRegression.PsychophysicalKernel
logits = data * weights(1:end-1) + weights(end);
log_bernoulli = responses(:) .* logits(:) - log(1 + exp(logits(:)));
LL = sum(log_bernoulli);
end

