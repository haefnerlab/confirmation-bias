function results = runSamplingModel(data, params)
%RUNSAMPLINGMODEL run the D->x->e sampling model.
%
% results = RUNSAMPLINGMODEL(data, params) gets model 'decisions' for data
% 'data' ([trials x frames] observed values of e). 'params' controls the
% behavior of the model. See Model.newSamplingParams for a descriptions of
% parameter options.
%
% The return value 'results' is a struct with the following fields:
%
%   results.params  - a copy of the 'params' argument
%   results.choices - [trials x 1] array of {-1, +1} values
%   results.x       - [trials x (frames*samples*batch+1)] sampled values of
%                     x (including extra 'initial' sample)
%   results.walk    - [trials x samples*frames+1] posterior log odds of
%                     D=+1/D=-1 (where walk(1) is the prior, hence size 
%                     frames+1)

if nargin < 2, params = Model.newSamplingParams(); end

%% Initialize return values

[trials, frames] = size(data);

sig_e = sqrt(params.var_e);
sig_x = sqrt(params.var_x);
p_match = params.p_match;
prior_D = params.prior_D;
gamma = params.gamma;
samples = params.samples;
batch = params.batch;

results = struct(...
    'params', params, ...
    'choices', zeros(trials, 1), ...
    'x', zeros(trials, frames*samples*batch+1), ...
    'walk', zeros(trials, frames*samples+1));

% Anonymous function to create (scaled) mixture-of-gaussians p(x|D)
p_x_D = @(D, pD) [+D sig_x p_match*pD; -D sig_x (1-p_match)*pD];
p_x_Dp = p_x_D(+1, 1); % 'p' for plus
p_x_Dm = p_x_D(-1, 1); % 'm' for minus

%% Run sampler in parallel
    function x = sample_x(e, pDp)
        mog_prior = [p_x_D(+1, pDp); p_x_D(-1, 1-pDp)];
        likelihood = [e, sig_e, 1];
        x = mogsample(mogprod(mog_prior, likelihood));
    end

% TODO - use parfor
for i=1:trials
    if mod(i, 10) == 0, disp(['trial ' num2str(i)]); end
    
    log_post_D = [log(prior_D) log(1-prior_D)];
    results.walk(i, 1) = log_post_D(1) - log_post_D(2);

    % draw initial samples using ancestral sampling
    D = sign(prior_D - rand); % this gives +1 with probability prior_D
    results.x(i, 1) = mogsample(p_x_D(D, 1));

    % Loop over evidence frames. Each frame, we apply 'samples' updates to
    % p(D), once for each 'batch' of samples from x.
    s_idx = 2;
    w_idx = 2;
    for j=1:frames
        e = data(i, j);
        for s=1:samples
            like_e_D = [0 0]; % sum of [p(e|D=+1) p(e|D=-1)] estimates
            % post_D contains [p(D=+1|e_1,...,e_j-1), p(D=-1|e_1,...,e_j-1)]
            post_D = exp(log_post_D - max(log_post_D)); post_D = post_D / sum(post_D);
            for b=1:batch
                % sample x using best posterior estimate of D
                results.x(i, s_idx) = sample_x(e, post_D(1));
                % accumulate estimate of p(e|D) using this sample of x
                like_e_D_cur = [mogpdf(results.x(i, s_idx), p_x_Dp), ...
                                mogpdf(results.x(i, s_idx), p_x_Dm)];
                like_e_D = like_e_D + like_e_D_cur / sum(like_e_D_cur);
                % increment sample index
                s_idx = s_idx + 1;
            end
            % update log_post_D with log of p(e|D), subtracting out the
            % previous log[p(D|e)]
            log_post_D = (1 - gamma / samples) * log_post_D + ...
                log(like_e_D / sum(like_e_D)) / samples;
            % record the posterior log odds in results.walk
            results.walk(i, w_idx) = log_post_D(1) - log_post_D(2);
            % increment walk index
            w_idx = w_idx + 1;
        end
    end
    
    results.choices(i) = sign(results.walk(i, end));
end

end

%% Helper functions for working with mixtures of gaussians (MoG)
% Note on format: a MoG (in 1d only) is specified by a mean, variance, and
% weight at each mode. A distribution with N modes is represented by an Nx3
% matrix:
%
%   mog = [mu_1 sigma_1 pi_1; ...; mu_n sigma_n, pi_n]

function x = mogsample(mog)
% first choose a mode
mode = find(rand <= cumsum(mog(:,3)), 1);
% sample from chosen mode
x = mog(mode, 1) + randn * mog(mode, 2);
end

function prod = mogprod(d1, d2)
modes1 = size(d1, 1);
modes2 = size(d2, 1);
modes_out = modes1 * modes2;

prod = zeros(modes_out, 3);
k = 1;
for i=1:modes1
    for j=1:modes2
        mu_i = d1(i,1); var_i = d1(i,2)^2; pi_i = d1(i,3);
        mu_j = d2(j,1); var_j = d2(j,2)^2; pi_j = d2(j,3);
        mu_k = (var_i * mu_j + var_j * mu_i) / (var_i + var_j);
        var_k = (var_i * var_j) / (var_i + var_j);
        pi_k = pi_i * pi_j * normpdf(mu_i, mu_j, sqrt(var_i + var_j));
        prod(k, :) = [mu_k, sqrt(var_k), pi_k];
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