function K = marginalJointKernel(XN, XM, theta)
%FITTING.MARGINALJOINTKERNEL Gaussian process kernel intended to be used to fit a log likelihood or
%log posterior over parameters. The idea is to separately estimate a GP over the (log) marginals,
%plus a GP over the joint. The latter is then expected to be simpler than if trying to just fit the
%joint directly.
%
%Parameters are [log(scale1) log(sigma1) ... log(scaleK) log(sigmaK) log(sigmaJoint)]

scale = exp(theta(1:2:end-1));
margStd = exp(theta(2:2:end-1));
jointStd = exp(theta(end));

K = zeros(size(XN,1), size(XM,1));

% Begin with sum of marginal covariances
for i=1:size(XN,2)
    K = K + margStd(i)^2 * exp(-(XN(:,i)-XM(:,i)').^2/(2*scale(i)^2));
end

% Add joint term
K = K + jointStd^2 * exp(-pdist2(XN, XM, 'seuclidean', scale).^2/2);
end