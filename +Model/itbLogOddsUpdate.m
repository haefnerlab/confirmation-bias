function [lpo] = itbLogOddsUpdate(params, e, lpo)
%MODEL.ISLOGODDSUPDATE compute update to log odds of C using importance sampling model.

noise = params.noise;
gamma = params.gamma;
bound = params.bound;

crossed = abs(lpo) >= bound;

% Only update for trials where bound hasn't been crossed yet
lpo(~crossed) = lpo(~crossed)*(1-gamma) + Model.logLikelihoodOdds(params, e(~crossed));

% Add zero-mean gaussian noise to lpo.
lpo = lpo + randn(size(lpo)) * noise;

% 'Stick' to bound after update
lpo = min(max(-bound, lpo), +bound);

end