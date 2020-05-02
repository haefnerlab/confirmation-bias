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
    val = Fitting.getParamsFields(params, field);
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

% simCountK is the counter, per trial, of how many simulations have been run for each trial since
% the last match.
simCountK = zeros(nTotalTrials, 1);

% Only keep re-running simulations on trials that haven't yet been matched 'repeats' times
matched = zeros(nTotalTrials, 1);

% At each iteration, the upper bound is based on L for all 'matched' trials plus what it would be
% for the 'unmatched' trials if they had all completed on the current iteration.
L_hat_upper_bound = zeros(1, repeats);

% Loop forever until no unmatched trials remain (but not an infinite loop b/c of upper/lower bound
% check)
while ~all(matched >= repeats)
    unmatched = matched < repeats;
    for iSet=length(params):-1:1
        if ~any(unmatched(setIdx(:,iSet))), continue; end
        
        % Run the model only on unmatched trials within this set
        sim_results = Model.runVectorized(params(iSet), signals(unmatched & setIdx(:,iSet), :) / params(iSet).signal_scale);
        
        % Count this simulation
        simCountK(unmatched & setIdx(:,iSet)) = simCountK(unmatched & setIdx(:,iSet)) + 1;
        
        % For all that now match... record 'k' and do book-keeping on what remains
        matching_subset = sim_results.choices == choices(unmatched & setIdx(:,iSet));
        if any(matching_subset)
            % Work out the indices out of 1:nTotalTrials that we just matched
            sim_idxs = find(unmatched & setIdx(:,iSet));
            match_idxs = sim_idxs(matching_subset);
            % Add 1 to 'matched' on matching trials
            matched(match_idxs) = matched(match_idxs) + 1;
            % Store and reset the 'k' for those trials
            for tr=match_idxs'
                next_nan = find(isnan(nSamplesToMatch(tr,:)),1);
                nSamplesToMatch(tr,next_nan) = simCountK(tr);
            end
            simCountK(match_idxs) = 0;
        end
    end
        
    % Enforce upper-bound on each repeat, i.e. on each column of 'nSamplesToMatch'. Upper-bound is
    % log likelihood estimate *if* all remaining simulations were matches.
    effSamples = nSamplesToMatch;
    effSimCount = simCountK;
    for r=1:repeats
        unmatched = isnan(nSamplesToMatch(:,r));
        % Skip columns that are done (note: includes those hit the lower bound already)
        if ~any(unmatched), continue; end
        
        % 'Effective' samples per trial set to hypothetical value of matching on the next iteration
        effSamples(unmatched, r) = effSimCount(unmatched)+1;

        % Maintaining the hypothetical, # sims per trial would drop to 1 for all remaining repeats
        effSimCount(unmatched) = 1;
    
        % See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
        % first derivative of log(gamma(x))
        L_hat_upper_bound(r) = sum(psi(0,1) - psi(0,effSamples(:,r)), 1);
        
        % If lower-bound has been crossed, set unmatched 'nSamplesToMatch' for repeat r to the
        % 'effective' samples from the upper-bound calculation.
        if L_hat_upper_bound(r) < lower_bound
            nSamplesToMatch(unmatched, r) = effSamples(unmatched, r);
            % Count these trials as matched.
            matched(unmatched) = matched(unmatched) + 1;
            % Reset 'simCountK' for these bound-crossed trials, again 'as if' it matches on the next
            % iteration.
            simCountK(unmatched) = 0;
        end
    end
end

% Sanity-check that all repeats have a value (even if set by crossing the bound)
assert(all(sum(~isnan(nSamplesToMatch), 2) == repeats) && all(nSamplesToMatch(:) > 0));
assert(all(matched == repeats));
    
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