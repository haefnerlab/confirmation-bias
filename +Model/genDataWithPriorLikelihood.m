function [data, var_e] = genDataWithPriorLikelihood(trials, frames, prior, likelihood)
%GENDATAWITHPRIORLIKELIHOOD generates a set of 'trials' (a 1xframes vector of real numbers),
%all with correct choice +1.

% From inverting the AUROC of two gaussians, with mu1=+1, mu2=-1, and
% identical variances:
%    p_auc = normcdf((mu1-mu2)/sqrt(2*var))
%          = normcdf(2 / (sqrt(2) * stdev))
% -> norminv(p_auc) = 2/sqrt(2) * 1/stdev
% -> stdev = sqrt(2) / norminv(p_auc)
stdev_e = sqrt(2) / norminv(likelihood);
var_e = stdev_e^2;

% generate the 'center' of each frame according to 'prior'; with
% probability 'prior' it is +1 and with probability '1-prior' it is -1.
centers = sign(prior - rand(trials, frames));

% draw signal from around the center with stdev calculated above.
data = normrnd(centers, stdev_e);
end