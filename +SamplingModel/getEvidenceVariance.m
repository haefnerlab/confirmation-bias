function var_e = getEvidenceVariance(p_sensory)
%GETEVIDENCEVARIANCE get variance of likelihood of P(e|x) given the desired 'sensory information'.
%
% AUROC for normal distributions is 
% 
% p = normcdf((mu_1 - mu_2) / sqrt(var_1 + var_2))
%
% Here, we solve for var_e=var_1=var_2 given that mu_1-mu_2=2. Note that norminv is the inverse of
% normcdf.
var_e = 2 ./ norminv(p_sensory).^2;
end