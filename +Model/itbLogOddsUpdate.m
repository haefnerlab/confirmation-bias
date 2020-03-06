function [lpo] = itbLogOddsUpdate(params, e, lpo)
%MODEL.ISLOGODDSUPDATE compute update to log odds of C using importance sampling model.

gamma = params.gamma;
bound = params.bound;

crossed = abs(lpo) >= bound;

% Only update for trials where bound hasn't been crossed yet
lpo(~crossed) = lpo(~crossed)*(1-gamma) + Model.logLikelihoodOdds(params, e(~crossed));

% 'Stick' to bound after upate
lpo = min(max(-bound, lpo), +bound);

end