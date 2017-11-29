function results = runSamplingModel(params)
%RUNSAMPLINGMODEL run the C->x->e sampling model.
%
% results = RUNSAMPLINGMODEL(params) gets model 'decisions' for an
% experiment with the given parameters. 'params' controls the behavior of
% the model and statistics of the data. See Model.newSamplingParams for a
% descriptions of parameter options.
%
% The return value 'results' is a struct with the following fields:
%
%   results.params  - a copy of the 'params' argument
%   results.choices - [trials x 1] array of {-1, +1} values
%   results.x       - [trials x (frames*samples*batch+1)] sampled values of
%                     x (including extra 'initial' sample)
%   results.w       - [trials x (frames*samples*batch+1)] weight of
%                     corresponding sample.
%   results.walk    - [trials x samples*frames+1] posterior log odds of
%                     C=+1/C=-1 (where walk(1) is the prior, hence size 
%                     frames+1)

if nargin < 1, params = Model.newModelParams(); end

data = Model.genDataWithParams(params);

%% Initialize return values

[trials, frames] = size(data);

sig_e = sqrt(params.var_e);
sig_x = sqrt(params.var_x);
p_match = params.p_match;
prior_C = params.prior_C;
samples = params.samples;
batch = params.batch;
gamma = params.gamma;

results = struct(...
    'params', params, ...
    'choices', zeros(trials, 1), ...
    'x', zeros(trials, frames*samples*batch+1), ...
    'w', zeros(trials, frames*samples*batch+1), ...
    'walk', zeros(trials, frames*samples+1));

% Create mixture-of-gaussian distributions
p = p_match * prior_C + (1 - p_match) * (1 - prior_C);
prior_x = mog.create([+1 -1], [sig_x sig_x], [p 1-p]);
p_x_Cp = mog.create([+1 -1], [sig_x sig_x], [p_match 1-p_match]); % 'p' for plus
p_x_Cm = mog.create([-1 +1], [sig_x sig_x], [p_match 1-p_match]); % 'm' for minus

%% Run sampler in parallel

% TODO - use parfor
for i=1:trials
    if mod(i, 10) == 0, disp(['trial ' num2str(i)]); end
    
    log_post_odds = log(prior_C) - log(1-prior_C);
    results.walk(i, 1) = log_post_odds;

    % draw initial samples using ancestral sampling
    results.x(i, 1) = mog.sample(prior_x);

    % Loop over evidence frames. Each frame, we apply 'samples' updates to
    % p(C), once for each 'batch' of samples from x.
    s_idx = 2;
    w_idx = 2;
    for j=1:frames
        e = data(i, j);
        for s=1:samples
            % post_C is current estimate of p(C=+1|e_1,...,e_j-1)
            post_C = exp(log_post_odds) / (1 + exp(log_post_odds));
            % Construct Q distribution to draw samples from
            likelihood = mog.create(e, sig_e, 1);
            prior = mog.create([+1 -1], [sig_x sig_x], [post_C 1-post_C]);
            Q = mog.prod(likelihood, prior);
            % sample x from Q 'batch' times and record result
            xs = mog.sample(Q, batch);
            results.x(i, s_idx:s_idx+batch-1) = xs;
            % Update accumulated evidence using importance-sampling weights
            ws = 1 ./ mog.pdf(xs, prior, params.importance_norm);
            results.w(i, s_idx:s_idx+batch-1) = ws;
            update = log(dot(mog.pdf(xs, p_x_Cp), ws)) - ...
                log(dot(mog.pdf(xs, p_x_Cm), ws));
            log_post_odds = (1 - gamma / samples) * log_post_odds + update / samples;
            % record the posterior log odds in results.walk
            results.walk(i, w_idx) = log_post_odds;
            % increment walk and sample index
            w_idx = w_idx + 1;
            s_idx = s_idx + batch;
        end
    end
    
    results.choices(i) = sign(results.walk(i, end));
end

end