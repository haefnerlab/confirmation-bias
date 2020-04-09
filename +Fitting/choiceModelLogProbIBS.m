function [log_post, log_like, est_variance, lower_bound, nSamplesToMatch] = choiceModelLogProbIBS(params, prior_info, signals, choices, lower_bound, repeats)
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

if nargin < 6 || isempty(repeats), repeats = 1; end

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
nSamplesToMatch = nan(nTotalTrials, repeats);

% Only keep re-running simulations on trials that haven't yet been matched 'repeats' times
matched = zeros(nTotalTrials, 1);

% At each iteration, the upper bound is based on L for all 'matched' trials plus what it would be
% for the 'unmatched' trials if they all completed on the next iteration.
L_hat_upper_bound = 0;

% Loop forever (up to maxK) until no unmatched trials remain
k = 1;
while ~all(matched >= repeats) && L_hat_upper_bound > lower_bound
    unmatched = matched < repeats;
    for iSet=length(params):-1:1
        if ~any(unmatched(setIdx(:,iSet))), continue; end
        
        % Run the model only on unmatched trials within this set
        sim_results = Model.runVectorized(params(iSet), signals(unmatched & setIdx(:,iSet), :) / params(iSet).signal_scale);
        
        % For all that now match... record 'k' and mark those trials as complete
        matching_subset = sim_results.choices == choices(unmatched & setIdx(:,iSet));
        if any(matching_subset)
            % Work out the indices out of 1:nTotalTrials that we just matched
            sim_idxs = find(unmatched & setIdx(:,iSet));
            match_idxs = sim_idxs(matching_subset);
            % Add 1 to 'matched' on matching trials
            matched(match_idxs) = matched(match_idxs) + 1;
            % Store the 'k' for this match number. If this is the 2nd or higher repeat for a given
            % trial, then 'k' is adjusted to number of simulations since the last success.
            last_match = nansum(nSamplesToMatch(match_idxs, :), 2);
            idx_update = sub2ind(size(nSamplesToMatch), match_idxs, matched(match_idxs));
            nSamplesToMatch(idx_update) = k-last_match;
        end
    end
        
    % Upper-bound on L-hat is the estimate of L under the assumption that all other trials will be
    % matched on the next step.
    effectiveSamples = nSamplesToMatch;
    effectiveSamples(matched == 0, 1) = k+1;
    
    % See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
    % first derivative of log(gamma(x))
    L_hat_upper_bound = sum(nanmean(psi(0,1) - psi(0,effectiveSamples), 2));

    % Advance...
    k = k+1;
end

if L_hat_upper_bound > lower_bound
    % If we didn't hit the lower bound, the following sanity-check should hold:
    assert(all(sum(~isnan(nSamplesToMatch), 2) == repeats));
else
    % We hit the lower bound.. use 'effectiveSamples' for the actual estimator
    nSamplesToMatch = effectiveSamples;
end
    
% See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
% first derivative of log(gamma(x))
log_lh_per_trial = nanmean(psi(0,1) - psi(0,nSamplesToMatch), 2);

% See eq. 16 in reference [1]; psi(1,x) is the trigamma function of x, i.e. the second
% derivative of log(gamma(x))
est_var_per_trial = nanmean(psi(1,1) - psi(1,nSamplesToMatch), 2) ./ repeats;
est_variance = sum(est_var_per_trial);

rng(randstate);

%% Combine prior and likelihood

% Unbiased estimate of the log posterior...
log_like = sum(log_lh_per_trial);
log_post = log_prior + log_like;
end