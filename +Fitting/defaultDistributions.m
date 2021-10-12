function distributions = defaultDistributions(fields, is_map)
distributions = struct(...
    'prior_C', betavar(2, 2, 20), ...
    'lapse', betavar(1, 10, 50), ...
    'gamma', uniform(0, 1, .1), ...
    'neggamma', uniform(-1, 1, .1), ...
    'signal_scale', gammavar(1, 20, 20), ...
    'step_size', gammavar(2, .05/2, 100), ...
    'var_x', gammavar(1, 1/10, 100), ...
    'noise', gammavar(1, 1/4, 100), ...
    'temperature', gammavar(1, 4, 100), ...gammavar(3/2, 1/3, 100), ...
    'bound', gammavar(2, 3, 60), ...
    'updates', integers(5, 10), ...
    'samples', integers(5, 20));

if nargin >= 1
    defaultFields = fieldnames(distributions);

    % Handle case where a field is prefixed with 'log_'
    islog = cellfun(@(f) startsWith(f, 'log_'), fields);
    fields = cellfun(@(f) strrep(f, 'log_', ''), fields, 'UniformOutput', false);

    % Handle case where field has an index, e.g. suffixed with '_1' or '_2'
    hasindex = cellfun(@(f) ~isempty(regexpi(f, '\w+\d*_\d+$')), fields);
    indexes = cellfun(@(f) regexprep(f, '\w+\d*(_\d+)', '$1'), fields, 'uniformoutput', false);
    indexes(~hasindex) = {''};
    fields = cellfun(@(f) regexprep(f, '_\d+$', ''), fields, 'UniformOutput', false);

    % Handle case where a field is suffixed with a number like  'lapse1' or 'lapse2' (e.g. multiple lapses)
    hassuffix = cellfun(@(f) ~isempty(regexpi(f, '\w+\d+$')), fields);
    suffixes = cellfun(@(f) regexprep(f, '\w+(\d+)', '$1'), fields, 'UniformOutput', false);
    suffixes(~hassuffix) = {''};
    fields = cellfun(@(f) regexprep(f, '\d+$', ''), fields, 'UniformOutput', false);

    % Warn about unrecognized fields and remove them
    unrecognized = setdiff(fields, defaultFields);
    if ~isempty(unrecognized)
        warning('The following fields are not recognized and have no default distributions: %s', strjoin(unrecognized, ', '));
        fields = setdiff(fields, unrecognized);
    end

    % Output struct will be built up one field at a time in the loop over fields below
    out_distribs = struct();

    % Post-process fields of the distributions struct with log transformations and name suffixes
    if nargin < 2, is_map = true; end
    for iF=1:length(fields)
        distrib = distributions.(fields{iF});
        full_name = fields{iF};
        if islog(iF)
            full_name = ['log_' full_name];
            olddistrib = distrib;

            % Adjust proposal random and bounds
            distrib.proprnd = @(logx, conc) log(olddistrib.proprnd(exp(logx), conc));
            distrib.priorrnd = @(sz) log(olddistrib.priorrnd(sz));
            distrib.lb  = log(olddistrib.lb);
            distrib.ub  = log(olddistrib.ub);
            distrib.plb = log(olddistrib.plb);
            distrib.pub = log(olddistrib.pub);

            % Adjust densities... Note that for true Bayesian inference, this must include a
            % change-of-variables adjustment |dx/dlogx|. However, for MAP inference we don't
            % want to change the location of the maximum through this adjustment. By default, assume
            % the MAP version. Let y=log(x). The change-of-variables adjustment is such that
            % P(y)dy=P(x)dx, so P(y)=P(e^y)dx/dy, or dx/dy=P(e^y) e^y
            if is_map
                % Simply transform logx --> x with no further density adjustments
                distrib.priorpdf = @(logx) olddistrib.priorpdf(exp(logx));
                distrib.proppdf = @(logx1,logx2,conc) olddistrib.proppdf(exp(logx1), exp(logx2), conc);
                distrib.logpriorpdf = @(logx) olddistrib.logpriorpdf(exp(logx));
                if isfield(distrib, 'logproppdf')
                    distrib.logproppdf = @(logx1,logx2,conc) olddistrib.logproppdf(exp(logx1), exp(logx2), conc);
                end
                
                % Having stated the following warning message, our strategy is to just ignore the
                % problem and use the "correct" change-of-variables adjusted iCDF here, even though
                % this is "wrong" from a MAP standpoint.
                warning('MAP behavior of iCDF on log-transformed var %s is not supported', upper(full_name));
                distrib.priorinv = @(u) log(olddistrib.priorinv(u));
            else
                distrib.priorpdf = @(logx) olddistrib.priorpdf(exp(logx)).*exp(logx);
                distrib.proppdf = @(logx1,logx2,conc) olddistrib.proppdf(exp(logx1),exp(logx2),conc).*exp(logx1);
                distrib.logpriorpdf = @(logx) olddistrib.logpriorpdf(exp(logx)) + logx;
                if isfield(distrib, 'logproppdf')
                    distrib.logproppdf = @(logx1,logx2,conc) olddistrib.logproppdf(exp(logx1), exp(logx2), conc) + logx1;
                end
                % Nothing fancy actually needed to do "correct" change of variables on the inverse
                % CDF function... just wrap it in a log() and done!
                distrib.priorinv = @(u) log(olddistrib.priorinv(u));
            end
        end

        full_name = [full_name suffixes{iF} indexes{iF}];

        out_distribs.(full_name) = distrib;
    end
    
    distributions = out_distribs;
end

end

function distrib = betavar(a, b, effective_pts)
% Create prior and proposal distributions for a beta-distributed variable with parameters [a,b].
% Proposal p(x1|x2) is based on treating the current sample as 'effective_pts' number of
% observations and sampling from the 'posterior'. Hence proposals concentrate around the current
% sample, with higher concentration for higher 'effective_pts'.
distrib = struct(...
    'priorpdf', @(x) betapdf(betaclip(x), a, b), ...
    'priorrnd', @(sz) betarnd(a, b, sz), ...
    'priorinv', @(u) betainv(u, a, b), ...
    'logpriorpdf', @(x) logbetapdf(betaclip(x), a, b), ...
    'proprnd', @(x,conc) betarnd(conc*(a + effective_pts * betaclip(x)), conc*(b + effective_pts * (1-betaclip(x)))), ...
    'proppdf', @(x1, x2, conc) betapdf(betaclip(x1), conc*(a + effective_pts * betaclip(x2)), conc*(b + effective_pts * (1-betaclip(x2)))), ...
    'logproppdf', @(x1, x2, conc) logbetapdf(betaclip(x1), conc*(a + effective_pts * betaclip(x2)), conc*(b + effective_pts * (1-betaclip(x2)))), ...
    'lb', betainv(.00001, a, b), 'plb', betainv(.01, a, b), 'ub', betainv(.99999, a, b), 'pub', betainv(.99, a, b));
end

function distrib = gammavar(k, t, prop_shape)
% Create prior and proposal distributions for a gamma-distributed variable with parameters [k,t].
% Proposal p(x1|x2) draws from a separate gamma distribution whose mean is the current parameter and
% whose shape is given by prop_shape. Higher prop_shape concentrates proposals around the current
% sample. 'k' is shape parameter: 1 makes it exponential, higher makes it more gaussian-like. Mean
% of the distribution is 'k*t'.
distrib = struct(...
    'priorpdf', @(x) gampdf(gammaclip(x), k, t), ...
    'priorrnd', @(sz) gamrnd(k, t, sz), ...
    'priorinv', @(u) gaminv(u, k, t), ...
    'logpriorpdf', @(x) loggampdf(gammaclip(x), k, t), ...
    'proprnd', @(x, conc) gamrnd(prop_shape*conc, gammaclip(x)/prop_shape/conc), ...
    'proppdf', @(x1, x2, conc) gampdf(gammaclip(x1), prop_shape*conc, gammaclip(x2)/prop_shape/conc), ...
    'logproppdf', @(x1, x2, conc) loggampdf(gammaclip(x1), prop_shape*conc, gammaclip(x2)/prop_shape/conc), ...
    'lb', gaminv(.00001, k, t), 'plb', gaminv(.01, k, t), 'ub', gaminv(.99999, k, t), 'pub', gaminv(.99, k, t));
end

function distrib =  integers(halfrange, tau)
% Create prior and proposal distribution for an integer-valued variable in {1, ..., inf}. Prior mass
% is proportional to an exponential distribution with mean tau. Proposals are generated uniformly
% between the current value plus or minus 'halfrange', clipped below at 1.
distrib = struct(...
    'priorpdf', @(i) exp(-i/tau), ...
    'priorrnd', @(sz) ceil(exprnd(tau, sz)), ...
    'priorinv', @(u) ceil(-tau*log(1-u)), ...
    'logpriorpdf', @(i) -i/tau, ...
    'proprnd', @(i,conc) randi([max(1, i-max(1,halfrange/sqrt(conc))), i+max(1,halfrange/sqrt(conc))]), ...
    'proppdf', @(i1, i2,conc) double(i1 >= max(1, i2-max(1,halfrange/sqrt(conc))) & i1 <= i2+max(1,halfrange/sqrt(conc))) / (i2+max(1,halfrange/sqrt(conc))-max(1,i2-max(1,halfrange/sqrt(conc)))+1), ...
    'lb', 1, 'plb', 1, 'ub', ceil(-tau*log(.00001)), 'pub', ceil(-tau*log(.01)));
end

function distrib = uniform(a, b, range)
% Create prior and proposal distribution for an real-valued variable in [a, b]. Proposals are
% generated uniformly between the current value plus or minus 'range'/2, clipped.
distrib = struct(...
    'priorpdf', @(x) ones(size(x))./(b-a), ...
    'priorrnd', @(sz) a+rand(sz)*(b-a), ...
    'priorinv', @(u) a*(1-u)+b*u, ...
    'logpriorpdf', @(x) -ones(size(x))*log(b-a), ...
    'proprnd', @(x,conc) max(a,x-range/2/sqrt(conc))+(min(b,x+range/2/sqrt(conc))-max(a,x-range/2/sqrt(conc)))*rand, ...
    'proppdf', @(x1, x2, conc) double(x1>=max(a,x2-range/2/sqrt(conc)) & x1<min(b,x2+range/2/sqrt(conc)))./(min(b,x2+range/2/sqrt(conc)) - max(a,x2-range/2/sqrt(conc))), ...
    'lb', a, 'plb', a+.01*(b-a), 'ub', b, 'pub', a+.99*(b-a));
end

function lp = logbetapdf(x, a, b)
lp = (a - 1) * log(x) + (b - 1) * log(1 - x) - betaln(a, b);
end

function lp = loggampdf(x, k, theta)
lp = (k - 1) * log(x) - x / theta - gammaln(k) - k * log(theta);
end

function val = betaclip(val)
val = min(max(val, 1e-9), 1-1e-9);
end
function val = gammaclip(val)
val = max(val, 1e-9);
end