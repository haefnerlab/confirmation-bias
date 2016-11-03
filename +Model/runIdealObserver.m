function results = runIdealObserver(data, params)
%RUNIDEALOBSERVER run an ideal observer on generative model specified by
%params.
%
% results = RUNIDEALOBSERVER(data, params) gets ideal 'decisions' and
% running posterior for the given data. Params are the same as the sampling
% params, but not all are used. The following are used:
%
%   params.var_e   - variance of gaussian p(e|D)
%   params.p_match - weight of modes, i.e. p(e|D) =
%                    p_match*N(D,var_e)+(1-p_match)*N(-D,var_e)
%   params.prior_D - prior probability of D=+1
%
% The return value 'results' is a struct with the following fields:
%
%   results.params  - a copy of the 'params' argument
%   results.choices - [trials x 1] array of {-1, +1} values
%   results.walk    - [trials x frames+1] posterior log odds of D=+1/D=-1
%                     (where walk(1) is the prior, hence size frames+1)

results = struct();
results.params = params;

% pDp is a mixture of gaussians, the probability that D is +1, and likewise
% pDm is for D=-1
pDp = [+1 sqrt(params.var_e) params.p_match;
       -1 sqrt(params.var_e) 1-params.p_match];
pDm = [-1 sqrt(params.var_e) params.p_match;
       +1 sqrt(params.var_e) 1-params.p_match];

% Compute log likelihood of each data point
log_odds = log(arrayfun(@(e) mogpdf(e, pDp), data)) ...
         - log(arrayfun(@(e) mogpdf(e, pDm), data));

% Walk is the log posterior odds. Log posterior is the cumuluative sum of
% the log prior and log likelihoods.
log_prior = log(params.prior_D) - log(1 - params.prior_D);
results.walk = cumsum([log_prior*ones(size(data, 1), 1) log_odds], 2);
results.choices = sign(results.walk(:, end));
end

function l = mogpdf(x, mog)
modes = size(mog, 1);
l_each_mode = normpdf(x * ones(modes, 1), mog(:,1), mog(:,2));
l = dot(l_each_mode, mog(:,3));
end