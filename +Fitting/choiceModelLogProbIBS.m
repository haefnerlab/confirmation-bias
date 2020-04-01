function [log_post, log_like, est_variance, lower_bound] = choiceModelLogProbIBS(params, prior_info, signals, choices, lower_bound)
% FITTING.choiceModelLogProbIBS estimate log posterior probability of given params struct. Prior is
% based on prior_info struct (see @Fitting.defaultDistributions). If prior_info is empty, prior is
% not counted and only the log likelihood is computed. The log likelihood is stochastic but unbiased,
% computed using the Inverse Binomial Sampling method of [1], clipped at a maximum of maxK
% iterations.
%
% Optionally, params may be a struct array of params, and both 'signals' and 'choices' may be cell
% arrays of the same size. The likelihood is then the product of likelihoods for each 'set' of data.
%
% [1] van Opheusden, B., Acerbi, L., & Ma, W. J. (2020). Unbiased and Efficient Log-Likelihood
%   Estimation with Inverse Binomial Sampling. http://arxiv.org/abs/2001.03985

if ~iscell(signals)
    signals = {signals};
    choices = {choices};
end

nSets = length(params);
assert(nSets == length(signals));
assert(nSets == length(choices));

if all(cellfun(@isempty, signals)), return; end

nTotalTrials = sum(cellfun(@length, choices));
if nargin < 5 || isempty(lower_bound)
    % By default, the lower bound assumes a uniform distribution over the two choices each trial.
    lower_bound = -nTotalTrials * log(2);
end

% Sanity check that all choices are in [-1 +1]
for iSet=1:nSets
    if ~all(ismember(choices{iSet}, [-1 +1]))
        error('Choices in wrong space! Should be in [-1 +1] but seeing %s.', mat2str(unique(choices{iSet})));
    end
end

% Out of the indices in 1:nTotalTrials, they are split into 'sets', one per triplet of
% (param,signal,choice). 'setIdx' is the logical mask of which trials belong to which set.
setIdx = false(nTotalTrials, nSets);
iNextStart = 1;
for iSet=1:nSets
    setIdx(iNextStart:iNextStart+length(choices{iSet})-1, iSet) = true;
    iNextStart = iNextStart+length(choices{iSet});
end

% Stack data vertically now that all 'sets' information has been processed
signals = vertcat(signals{:});
choices = vertcat(choices{:});

%% Evaluate prior
log_prior = 0;
prior_fields = fieldnames(prior_info);
for iF=1:length(prior_fields)
    field = prior_fields{iF};
    val = Fitting.getParamsFields(params(1), field);
    if isfield(prior_info.(field), 'logpriorpdf')
        log_prior = log_prior + prior_info.(field).logpriorpdf(val);
    else
        log_prior = log_prior + log(prior_info.(field).priorpdf(val));
    end
end

%% Evaluate likelihood using Inverse Binomial Sampling

% Maximize pseudo-randomness inside the simulations
randstate = rng();
rng('shuffle');

% All we need for the IBS likelihood is the number of model simulations, each trial, until we
% match the choice on that trial
nSamplesToMatch = ones(nTotalTrials, 1);

% Only keep re-running simulations on trials that haven't yet been matched
matched = false(nTotalTrials, 1);

% At each iteration, the upper bound is based on L for all 'matched' trials plus what it would be
% for the 'unmatched' trials if they all completed on the next iteration.
L_hat_upper_bound = 0;

% Loop forever (up to maxK) until no unmatched trials remain
k = 1;
while ~all(matched) && L_hat_upper_bound > lower_bound
    for iSet=length(params):-1:1
        if all(matched(setIdx(:,iSet))), continue; end
        
        % Run the model only on unmatched trials within this set
        sim_results = Model.runVectorized(params(iSet), signals(~matched & setIdx(:,iSet), :));
        
        % For all that now match... record 'k' and mark those trials as complete
        matching_subset = sim_results.choices == choices(~matched & setIdx(:,iSet));
        if any(matching_subset)
            % Work out the indices out of 1:nTotalTrials that we just matched
            sim_idxs = find(~matched & setIdx(:,iSet));
            match_idxs = sim_idxs(matching_subset);
            % Mark those trials as completed
            matched(match_idxs) = true;
            % Store their 'k'
            nSamplesToMatch(match_idxs) = k;
        end
    end
        
    % Upper-bound on L-hat is the estimate of L under the assumption that all other trials will be
    % matched on the next step.
    nSamplesToMatch(~matched) = k+1;
    
    % See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
    % first derivative of log(gamma(x))
    L_hat_upper_bound = sum(psi(0,1) - psi(0,nSamplesToMatch));

    % Advance...
    k = k+1;
end
    
% See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
% first derivative of log(gamma(x))
log_lh_per_trial = psi(0,1) - psi(0,nSamplesToMatch);

% See eq. 16 in reference [1]; psi(1,x) is the trigamma function of x, i.e. the second
% derivative of log(gamma(x))
var_ll_per_trial = psi(1,1) - psi(1,nSamplesToMatch);

rng(randstate);

%% Combine prior and likelihood

% Unbiased estimate of the log posterior...
log_like = sum(log_lh_per_trial);
log_post = log_prior + log_like;

% Estimate of the variance in log_post
est_variance = sum(var_ll_per_trial);
end