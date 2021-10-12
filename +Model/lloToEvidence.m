function s = lloToEvidence(var_s, sensory_llos)
% Convert LLOs to model evidence scalars by solving the relation
%
%     log(normpdf(s,mu_1,sig_1)/normpdf(s,mu_2,sig_2))=llo
%
% In general, this would be a quadratic in s. Here, we know sig_1=sig_2, and mu_1=-mu_2=+1, which
% greatly simplifies things. The log likelihood odds, or 'llo' must be 'sensory only', i.e. each
% distribution is a single gaussian rather than a mixture of them.

s = sensory_llos * var_s / 4;

end