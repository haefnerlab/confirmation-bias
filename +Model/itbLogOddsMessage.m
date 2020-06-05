function [pr_lpo] = itbLogOddsMessage(params, e, lpo_edges, pr_lpo)
%MODEL.ISLOGODDSMESSAGE same external behavior as @Model.itbLogOddsUpdate, but rather than noise
%draws (randn), propagates a discrete distribution of pr(lpo). Input lpo_edges specifies the
%location of the edges of bins, where pr_lpo is [trials x 1 x bins] probability mass. Expects
%edges(1)=-inf and edges(end)=+inf so that the 1st and last bins hold mass of bound-crossings.

%Note: looping over trials is slow. Proper vectorization still to-do.

noise = params.noise;
gamma = params.gamma;
bound = params.bound;

% Precondition on lpo_edges is that the endpoints contain the integration bounds.
lower_bound_contained = lpo_edges(1) <= -bound && -bound <= lpo_edges(2);
upper_bound_contained = lpo_edges(end-1) <= +bound && +bound <= lpo_edges(end);
assert(isinf(bound) || (lower_bound_contained && upper_bound_contained), 'Edges do not contain bounds!');

% Enforce sanity: mass per trial must sum to 1
pr_lpo = pr_lpo ./ sum(pr_lpo, 3);

% Compute log likelihood odds given e, which is the amount of 'shift' to mass per trial
llo = Model.logLikelihoodOdds(params, e);

for tr=1:length(e)
    pmf = pr_lpo(tr, 1, :);
    
    % Hold on to out-of-bounds mass
    lo = pmf(1); hi = pmf(end);
    
    % Re-normalize mass in the 'middle' that is still being integrated (hence "int_")
    int_pmf = pmf(2:end-1);
    int_mass = sum(pmf(2:end-1));
    int_pmf = int_pmf(:)' / int_mass;
    
    if int_mass < 1e-9
        % All mass is in the bounds. Don't bother doing any more integration.
        pr_lpo(tr, 1, 1) = lo;
        pr_lpo(tr, 1, end) = hi;
        pr_lpo(tr, 1, 2:end-1) = 0;
        continue;
    end
    
    % Propagate forward, shifting by llo, scaling by (1-gamma), and adding noise
    [new_pmf, new_lo, new_hi] = Model.propagateNoisePMF(lpo_edges(2:end-1), int_pmf, llo(tr), 1-gamma, noise);
    
    % Accumulate total bound crossings and rescale 'center' part by initial amount of mass there
    pr_lpo(tr, 1, 1) = lo + new_lo * int_mass;
    pr_lpo(tr, 1, end) = hi + new_hi * int_mass;
    pr_lpo(tr, 1, 2:end-1) = new_pmf * int_mass;
end

end