function distributions = defaultDistributions(fields, is_map, allow_gamma_neg)
distributions = struct(...
    'prior_C', betavar(2, 2, 50), ...
    'lapse', betavar(1, 10, 50), ...
    'gamma', uniform(0, 1, .05), ...
    'sensor_noise', gammavar(1, 3, 10), ...
    'var_x', gammavar(1, 1/10, 10), ...
    'noise', gammavar(1, 1/4, 10), ...
    'temperature', gammavar(1, 8, 10), ...gammavar(3/2, 1/3, 10), ...
    'bound', gammavar(2, 3, 10), ...
    'updates', integers(5, 10), ...
    'samples', integers(5, 20));

if nargin >= 3 && allow_gamma_neg
    distributions.gamma = uniform(-1, 1, .05);
end

if nargin >= 1
    defaultFields = fieldnames(distributions);

    % Handle case where a field is prefixed with 'log_'
    islog = cellfun(@(f) startsWith(f, 'log_'), fields);
    fields = cellfun(@(f) strrep(f, 'log_', ''), fields, 'UniformOutput', false);

    % Handle case where a field is suffixed with a number like  '_1' or '_2' (e.g. multiple lapses)
    hassuffix = cellfun(@(f) ~isempty(regexpi(f, '_\d+$')), fields);
    suffixes = cellfun(@(f) regexprep(f, '.*(?=_\d+$)', ''), fields, 'UniformOutput', false);
    suffixes(~hassuffix) = {''};
    fields = cellfun(@(f) regexprep(f, '_\d+$', ''), fields, 'UniformOutput', false);

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
            distrib.proprnd = @(logx) log(olddistrib.proprnd(exp(logx)));
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
                distrib.proppdf = @(logx) olddistrib.proppdf(exp(logx));
                distrib.logpriorpdf = @(logx) olddistrib.logpriorpdf(exp(logx));
                if isfield(distrib, 'logproppdf')
                    distrib.logproppdf = @(logx) olddistrib.logproppdf(exp(logx));
                end
            else
                distrib.priorpdf = @(logx) olddistrib.priorpdf(exp(logx)).*exp(logx);
                distrib.proppdf = @(logx) olddistrib.proppdf(exp(logx)).*exp(logx);
                distrib.logpriorpdf = @(logx) olddistrib.logpriorpdf(exp(logx)) + logx;
                if isfield(distrib, 'logproppdf')
                    distrib.logproppdf = @(logx) olddistrib.logproppdf(exp(logx)) + logx;
                end
            end
        end

        if hassuffix(iF)
            full_name = [full_name suffixes{iF}];
        end

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
    'logpriorpdf', @(x) logbetapdf(betaclip(x), a, b), ...
    'proprnd', @(x) betarnd(a + effective_pts * betaclip(x), b + effective_pts * (1-betaclip(x))), ...
    'proppdf', @(x1, x2) betapdf(betaclip(x1), a + effective_pts * betaclip(x2), b + effective_pts * (1-betaclip(x2))), ...
    'logproppdf', @(x1, x2) logbetapdf(betaclip(x1), a + effective_pts * betaclip(x2), b + effective_pts * (1-betaclip(x2))), ...
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
    'logpriorpdf', @(x) loggampdf(gammaclip(x), k, t), ...
    'proprnd', @(x) gamrnd(prop_shape, gammaclip(x)/prop_shape), ...
    'proppdf', @(x1, x2) gampdf(gammaclip(x1), prop_shape, gammaclip(x2)/prop_shape), ...
    'logproppdf', @(x1, x2) loggampdf(gammaclip(x1), prop_shape, gammaclip(x2)/prop_shape), ...
    'lb', gaminv(.00001, k, t), 'plb', gaminv(.01, k, t), 'ub', gaminv(.99999, k, t), 'pub', gaminv(.99, k, t));
end

function distrib =  integers(halfrange, tau)
% Create prior and proposal distribution for an integer-valued variable in {1, ..., inf}. Prior mass
% is proportional to an exponential distribution with mean tau. Proposals are generated uniformly
% between the current value plus or minus 'halfrange', clipped below at 1.
distrib = struct(...
    'priorpdf', @(i) exp(-i/tau), ...
    'priorrnd', @(sz) ceil(exprnd(tau, sz)), ...
    'logpriorpdf', @(i) -i/tau, ...
    'proprnd', @(i) randi([max(1, i-halfrange), i+halfrange]), ...
    'proppdf', @(i1, i2) double(i1 >= max(1, i2-halfrange) && i1 <= i2+halfrange) / (i2+halfrange-max(1,i2-halfrange)+1), ...
    'lb', 1, 'plb', 1, 'ub', -tau*log(.00001), 'pub', -tau*log(.01));
end

function distrib = uniform(a, b, range)
% Create prior and proposal distribution for an real-valued variable in [a, b]. Proposals are
% generated uniformly between the current value plus or minus 'range'/2, clipped.
distrib = struct(...
    'priorpdf', @(x) ones(size(x))./(b-a), ...
    'priorrnd', @(sz) a+rand(sz)*(b-a), ...
    'logpriorpdf', @(x) -ones(size(x))*log(b-a), ...
    'proprnd', @(x) max(a,x-range/2)+(min(b,x+range/2)-max(a,x-range/2))*rand, ...
    'proppdf', @(x1, x2) 1./(min(b,x+range/2) - max(a,x-range/2)), ...
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