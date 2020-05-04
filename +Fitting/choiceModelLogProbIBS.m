function [log_post, log_like, est_variance, lower_bound, nSamplesToMatch] = choiceModelLogProbIBS(params, prior_info, signals, choices, lower_bound, repeats, bailoutK)
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
if nargin < 7 || isempty(bailoutK), bailoutK = inf; end

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

% As in @Fitting.choiceModelLogProb, we keep track of the average 'prob_choice' output by the model
% as well as the variance in this value across simulations. If K > bailoutK, we revert to using an
% estimator of log(p) based on a variance-corrected log(avgP)
avgP = zeros(nTotalTrials, 1);
varP_S = zeros(nTotalTrials, 1);
bailout_ll = -inf(nTotalTrials, repeats);
bailout_ll_var = inf(nTotalTrials, repeats);

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
        
        sim_idx = unmatched & setIdx(:,iSet);
        
        % Run the model only on unmatched trials within this set
        sim_results = Model.runVectorized(params(iSet), signals(sim_idx, :) / params(iSet).signal_scale);
        
        % Count this simulation
        simCountK(sim_idx) = simCountK(sim_idx) + 1;
        
        % Update mean and variance trackers of simulation probability. Both of these are 'running'
        % estimators for mean and variance of p per trial. See
        % https://www.johndcook.com/blog/standard_deviation/
        newP = sim_results.prob_choice;
        newAvgP = avgP(sim_idx)+(newP-avgP(sim_idx))./simCountK(sim_idx);
        varP_S(sim_idx) = varP_S(sim_idx) + (newP-avgP(sim_idx)).*(newP-newAvgP);
        avgP(sim_idx) = newAvgP;
        
        % For all that now match... record 'k' and do book-keeping on what remains
        matching_subset = sim_results.choices == choices(sim_idx);
        if any(matching_subset)
            % Work out the indices out of 1:nTotalTrials that we just matched
            sim_idxs = find(sim_idx);
            match_idxs = sim_idxs(matching_subset);
            % Add 1 to 'matched' on matching trials
            matched(match_idxs) = matched(match_idxs) + 1;
            % Store and reset the 'k' for those trials
            for tr=match_idxs'
                r = find(isnan(nSamplesToMatch(tr,:)),1);
                nSamplesToMatch(tr,r) = simCountK(tr);
            end
            simCountK(match_idxs) = 0;
            avgP(match_idxs) = 0;
            varP_S(match_idxs) = 0;
        end
    end
    
    % Check for bailout condition - revert to average-of-prob estimator
    if any(simCountK >= bailoutK)
        bailoutTr = find(simCountK >= bailoutK);
        for tr=bailoutTr'
            r = find(isnan(nSamplesToMatch(tr,:)),1);
            
            % Record a precise albeit biased estimate of the true ll for this trial. First, convert
            % from varP_S counter to an actual estimate of var(p) for this trial:
            varP = varP_S(tr) / (simCountK(tr)-1);
            % Next, use estimate of mean of p (muP) and its variance (varP) to construct an
            % estimate of log(p) and its variance. See stats.stackexchange.com/a/57766/234036
            muP = avgP(tr).*(choices(tr)==+1) + (1-avgP(tr)).*(choices(tr)==-1);
            bailout_ll(tr, r) = log(muP) - 1/2*varP/muP^2;
            bailout_ll_var(tr, r) = varP./muP^2;
            % Issue a warning if these estimates are likely to be off by a lot. See comments in
            % stats.stackexchange.com/a/57766/234036
            if muP/sqrt(varP) < 1.5
                warning('SNR of value bailout calculation is too low for confident estimates of LL or var(LL)!');
            elseif muP/sqrt(varP) < 2.5
                warning('SNR of value bailout calculation is too low for confident estimates of var(LL)!');
            end
            
            % Count this as a match after 'bailoutK' simulations (but this will not actually be used
            % for the LL estimate)
            nSamplesToMatch(tr,r) = bailoutK;
        end
        matched(bailoutTr) = matched(bailoutTr)+1;
        simCountK(bailoutTr) = 0;
        avgP(bailoutTr) = 0;
        varP_S(bailoutTr) = 0;
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

        % Maintaining the hypothetical, # sims per trial would drop to 0 after a match
        effSimCount(unmatched) = 0;
    
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
            avgP(unmatched) = 0;
            varP_S(unmatched) = 0;
        end
    end
end

% Sanity-check that all repeats have a value (even if set by crossing the bound)
assert(all(sum(~isnan(nSamplesToMatch), 2) == repeats) && all(nSamplesToMatch(:) > 0));
assert(all(matched == repeats));
    
% See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
% first derivative of log(gamma(x))
log_lh_estimates = psi(0,1) - psi(0,nSamplesToMatch);
log_lh_estimates(nSamplesToMatch >= bailoutK) = bailout_ll(nSamplesToMatch >= bailoutK);
log_lh_per_trial = nanmean(log_lh_estimates, 2);

% See eq. 16 in reference [1]; psi(1,x) is the trigamma function of x, i.e. the second
% derivative of log(gamma(x))
var_log_lh_estimates = psi(1,1) - psi(1,nSamplesToMatch);
var_log_lh_estimates(nSamplesToMatch >= bailoutK) = bailout_ll_var(nSamplesToMatch >= bailoutK);
est_var_per_trial = nanmean(var_log_lh_estimates, 2) ./ repeats;
est_variance = sum(est_var_per_trial);

rng(randstate);

%% Combine prior and likelihood

% Unbiased estimate of the log posterior...
log_like = sum(log_lh_per_trial);
log_post = log_prior + log_like;
end