function var_s = getEvidenceVariance(p_sensory)
%GETEVIDENCEVARIANCE get variance of likelihood of P(s|x) given the desired 'sensory information'.
%
% AUROC for normal distributions is 
% 
% p = normcdf((mu_1 - mu_2) / sqrt(var_1 + var_2))
%
% Here, we solve for var_s=var_1=var_2 given that mu_1-mu_2=2. Note that norminv is the inverse of
% normcdf.
var_s = 2 ./ norminv(p_sensory).^2;
end