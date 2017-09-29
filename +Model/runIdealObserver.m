function results = runIdealObserver(params)
%RUNIDEALOBSERVER run an ideal observer on generative model specified by
%params.
%
% results = RUNIDEALOBSERVER(data, params) gets ideal 'decisions' and
% running posterior for the given data. Params are the same as the sampling
% params, but not all are used. The following are used:
%
%   params.var_e   - variance of gaussian p(e|C)
%   params.p_match - weight of modes, i.e. p(e|C) =
%                    p_match*N(C,var_e)+(1-p_match)*N(-C,var_e)
%   params.prior_C - prior probability of C=+1
%
% The return value 'results' is a struct with the following fields:
%
%   results.params  - a copy of the 'params' argument
%   results.choices - [trials x 1] array of {-1, +1} values
%   results.walk    - [trials x frames+1] posterior log odds of C=+1/C=-1
%                     (where walk(1) is the prior, hence size frames+1)

data = Model.genDataWithParams(params);
results = struct();
results.params = params;

p_match = params.p_match;
stdev = sqrt(params.var_e);

% pCp is a mixture of gaussians, the probability that C is +1, and likewise
% pCm is for C=-1
pCp = mog.create([+1 -1], [stdev stdev], [p_match 1-p_match]);
pCm = mog.create([-1 +1], [stdev stdev], [p_match 1-p_match]);

% Compute log likelihood of each data point
log_odds = mog.logpdf(data, pCp) - mog.logpdf(data, pCm);

% Walk is the log posterior odds. Log posterior is the cumuluative sum of
% the log prior and log likelihoods.
log_prior = log(params.prior_C) - log(1 - params.prior_C);
results.walk = cumsum([log_prior*ones(size(data, 1), 1) log_odds], 2);
results.choices = sign(results.walk(:, end));
end