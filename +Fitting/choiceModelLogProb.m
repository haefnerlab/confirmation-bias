function [log_post, log_prior, log_lh_per_trial, est_variance] = choiceModelLogProb(params, prior_info, signals, choices, maxK)
% FITTING.CHOICEMODELLOGPROB estimate log posterior probability of given params struct. Prior is
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

if nargin < 5, maxK = 1e4; end

if ~iscell(signals)
    signals = {signals};
    choices = {choices};
end

assert(length(params) == length(signals));
assert(length(params) == length(choices));

log_lh_per_trial = cell(size(params));

if all(cellfun(@isempty, signals)), return; end

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

log_lh = 0;
for iSet=length(params):-1:1
    nTrials = length(choices{iSet});
    assert(all(ismember(choices{iSet}, [-1 +1])));

    % All we need for the IBS likelihood is the number of model simulations, each trial, until we
    % match the choice on that trial
    nSamplesToMatch = maxK * ones(nTrials, 1);

    % Only keep re-running simulations on trials that haven't yet been matched
    unmatched = 1:nTrials;
    
    % Loop forever (up to maxK) until no unmatched trials remain
    k = 1;
    while ~isempty(unmatched) && k < maxK
        % Run the model on unmatched trials
        sim_results = Model.runVectorized(params(iSet), signals{iSet}(unmatched, :));

        % For all that now match... record 'k' and remove those trials
        matching_subset = sim_results.choices == choices{iSet}(unmatched);
        nSamplesToMatch(unmatched(matching_subset)) = k;
        unmatched(matching_subset) = [];

        % Advance...
        k = k+1;
    end

    % See eq. 14 in reference [1]; psi(0,x) is Matlab's builtin digamma function of x, i.e. the
    % first derivative of log(gamma(x))
    log_lh_per_trial = psi(0,1) - psi(0,nSamplesToMatch);
    log_lh(iSet) = sum(log_lh_per_trial);
    
    % See eq. 16 in reference [1]; psi(1,x) is the trigamma function of x, i.e. the second
    % derivative of log(gamma(x))
    var_ll_per_trial{iSet} = psi(1,1) - psi(1,nSamplesToMatch);
end

rng(randstate);

%% Combine prior and likelihood

% Unbiased estimate of the log posterior...
log_post = log_prior + sum(log_lh);

% Estimate of the variance in log_post
all_ll_var = vertcat(var_ll_per_trial{:});
est_variance = sum(all_ll_var);
end