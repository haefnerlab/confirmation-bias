function [pr_lpo] = itbLogOddsMessage(params, e, lpo_edges, pr_lpo)
%MODEL.ISLOGODDSMESSAGE same external behavior as @Model.itbLogOddsUpdate, but rather than noise
%draws (randn), propagates a discrete distribution of pr(lpo). Input lpo_edges specifies the
%location of the edges of bins, where pr_lpo is [trials x bins] probability mass. Expects
%edges(1)=-inf and edges(end)=+inf so that the 1st and last bins hold mass of bound-crossings.

noise = params.noise;
gamma = params.gamma;
bound = params.bound;

% Precondition on lpo_edges is that the endpoints contain the integration bounds.
lower_bound_contained = lpo_edges(1) <= -bound && -bound <= lpo_edges(2);
upper_bound_contained = lpo_edges(end-1) <= +bound && +bound <= lpo_edges(end);
assert(isinf(bound) || (lower_bound_contained && upper_bound_contained), 'Edges do not contain bounds!');

% Enforce sanity: mass per trial must sum to 1
pr_lpo = pr_lpo ./ sum(pr_lpo, 2);

% Compute log likelihood odds given e, which is the amount of 'shift' to mass per trial
llo = Model.logLikelihoodOdds(params, e);

% Re-normalize mass in the 'middle' that is still being integrated (hence "int_")
int_pr_lpo = pr_lpo(:,2:end-1);
int_mass = sum(int_pr_lpo, 2);
int_pr_lpo = int_pr_lpo ./ int_mass;

% Only bother working with the trials with some probability mass that hasn't crossed bounds yet.
to_propagate = int_mass > 1e-9;
if ~any(to_propagate)
    return;
end

% Propagate forward, shifting by llo, scaling by (1-gamma), and adding noise
[new_pr_lpo, new_lo, new_hi] = Model.propagateNoisePMF(lpo_edges(2:end-1), int_pr_lpo(to_propagate,:), llo(to_propagate), 1-gamma, noise);

% Accumulate total bound crossings and rescale 'center' part by initial amount of mass there
pr_lpo(to_propagate, 1) = pr_lpo(to_propagate, 1) + new_lo .* int_mass(to_propagate);
pr_lpo(to_propagate, end) = pr_lpo(to_propagate, end) + new_hi .* int_mass(to_propagate);
pr_lpo(to_propagate, 2:end-1) = new_pr_lpo .* int_mass(to_propagate);

end