function results = runSamplingModel(data, params)
%RUNSAMPLINGMODEL run the D->x->e sampling model.
%
% results = RUNSAMPLINGMODEL(data, params) gets model 'decisions' for data
% 'data' ([trials x frames] observed values of e). 'params' controls the
% behavior of the model:
%
%   params.var_e   - variance of gaussian p(e|x)
%   params.var_x   - variance of gaussian(s) p(x|D)
%   params.p_x     - weight of x modes, i.e. p(x|D) =
%                    p_x*N(D,var_x)+(1-p_x)*N(-D,var_x)
%   params.prior_D - prior probability of D=+1
%   params.gamma   - value in [0,1]; how much the prior is 'subtracted
%                    out'. 0 = low variance, 1 = unbiased.
%   params.samples - num samples per piece of evidence
%
% The return value 'results' is a struct with the following fields:
%
%   results.params  - a copy of the 'params' argument
%   results.choices - [trials x 1] array of {-1, +1} values
%   results.x       - [trials x (frames*samples+1)] sampled values of x
%                     (including extra 'initial' sample)
%   results.D       - [trials x (frames*samples+1)] sampled values of D
%                     (including extra 'initial' sample)
%   results.walk    - [trials x frames+1] posterior log odds of D=+1/D=-1
%                     (where walk(1) is the prior, hence size frames+1)

%% Initialize return values
if nargin < 2, params = struct(); end
if ~isfield(params, 'var_e'), params.var_e = 0.1; end
if ~isfield(params, 'var_x'), params.var_x = 0.5; end
if ~isfield(params, 'p_x'), params.p_x = 1; end
if ~isfield(params, 'prior_D'), params.prior_D = 0.5; end
if ~isfield(params, 'gamma'), params.gamma = 0; end
if ~isfield(params, 'samples'), params.samples = 1; end

[trials, frames] = size(data);

var_e = params.var_e;
var_x = params.var_x;
p_x = params.p_x;
prior_D = params.prior_D;
gamma = params.gamma;
samples = params.samples;

results = struct(...
    'params', params, ...
    'choices', zeros(trials, 1), ...
    'x', zeros(trials, frames*samples+1), ...
    'D', zeros(trials, frames*samples+1), ...
    'walk', zeros(trials, frames+1));

% Anonymous function to create mixture-of-gaussians p(x|D)
p_x_D = @(D) [+D var_x p_x; -D var_x 1-p_x];
p_x_Dp = p_x_D(+1); % 'p' for plus
p_x_Dm = p_x_D(-1); % 'm' for minus

%% Run sampler in parallel

    function x = sample_x(e, D)
        mog_prior = p_x_D(D);
        likelihood = [e, var_e, 1];
        x = mogsample(mogprod(mog_prior, likelihood));
    end

    function D = sample_D(x, post_D)
        % pDp is "propability of D plus one", and pDm is "probability of D
        % minus one"
        pDp = post_D * mogpdf(x, p_x_Dp);
        pDm = (1 - post_D) * mogpdf(x, p_x_Dm);
        pDp = pDp / (pDp + pDm);
        D = sign(pDp - rand); % this gives +1 with probability pDp
    end

% TODO - use parfor
for i=1:trials
    if mod(i, 10) == 0, disp(['trial ' num2str(i)]); end
    
    log_post_D = [log(prior_D) log(1-prior_D)];
    results.walk(i, 1) = log_post_D(1) - log_post_D(2);

    % draw initial samples using ancestral sampling
    results.D(i, 1) = sign(prior_D - rand); % this gives +1 with probability prior_D
    results.x(i, 1) = mogsample(p_x_D(results.D(i, 1)));

    % loop over evidence frames
    s_idx = 2;
    for j=1:frames
        e = data(i, j);
        like_e_D = [0 0]; % sum of [p(e|D=+1) p(e|D=-1)] estimates
        % post_D contains [p(D=+1|e_1,...,e_j-1), p(D=-1|e_1,...,e_j-1)]
        post_D = exp(log_post_D); post_D = post_D / sum(post_D);
        for s=1:samples
            % sample x conditioned on the previous D
            results.x(i, s_idx) = sample_x(e, results.D(i, s_idx-1));
            % sample D conditioned on the just-sampled x
            results.D(i, s_idx) = sample_D(results.x(i, s_idx), post_D);
            % accumulate estimate of p(e|D) using current x
            like_e_D_cur = [mogpdf(results.x(i, s_idx), p_x_Dp), ...
                            mogpdf(results.x(i, s_idx), p_x_Dm)];
            like_e_D = like_e_D + like_e_D_cur / sum(like_e_D_cur);
            % increment sample index
            s_idx = s_idx + 1;
        end
        % update log_post_D with log of p(e|D), subtracting out the
        % previous log[p(D|e)]
        log_post_D = (1 - gamma) * log_post_D + log(like_e_D / sum(like_e_D));
        % record the posterior log odds in results.walk
        results.walk(i, j+1) = log_post_D(1) - log_post_D(2);
    end
    
    results.choices(i) = sign(results.walk(i, end));
end

end

%% Helper functions for working with mixtures of gaussians (MoG)
% Note on format: a MoG (in 1d only) is specified by a mean, variance, and
% weight at each mode. A distribution with N modes is represented by an Nx3
% matrix:
%
%   mog = [mu_1 var_1 pi_1; ...; mu_n var_n, pi_n]

function x = mogsample(mog)
% first choose a mode
mode = find(rand <= cumsum(mog(:,3)), 1);
% sample from chosen mode
x = normrnd(mog(mode, 1), mog(mode, 2));
end

function prod = mogprod(d1, d2)
modes1 = size(d1, 1);
modes2 = size(d2, 1);
modes_out = modes1 * modes2;

prod = zeros(modes_out, 3);
k = 1;
for i=1:modes1
    for j=1:modes2
        mu_i = d1(i,1); var_i = d1(i,2); pi_i = d1(i,3);
        mu_j = d2(j,1); var_j = d2(j,2); pi_j = d2(j,3);
        mu_k = (var_i * mu_j + var_j * mu_i) / (var_i + var_j);
        var_k = (var_i * var_j) / (var_i + var_j);
        pi_k = pi_i * pi_j * normpdf(mu_i, mu_j, sqrt(var_i + var_j));
        prod(k, :) = [mu_k, var_k, pi_k];
        k = k+1;
    end
end
prod(:,3) = prod(:,3) / sum(prod(:,3));
end

function l = mogpdf(x, mog)
modes = size(mog, 1);
l_each_mode = normpdf(x * ones(modes, 1), mog(:,1), mog(:,2));
l = dot(l_each_mode, mog(:,3));
end