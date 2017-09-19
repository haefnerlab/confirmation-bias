function plotISBiasVariance(pri, like, sampleses, n_gridpts, n_repeats, use_precomputed)

if nargin < 4, use_precomputed = false; end

%% Load or compute results
savedir = fullfile('+Model', 'saved results');

log_prior_c = linspace(-2, 2, n_gridpts);
e_values = linspace(-2, 2, n_gridpts);
gamma_values = linspace(0, 1, n_gridpts);
[ee, pp, gg] = meshgrid(e_values, log_prior_c, gamma_values);

bias2_results = cell(size(sampleses));
variance_results = cell(size(sampleses));
marginal_bias2 = cell(size(sampleses));
marginal_variance = cell(size(sampleses));

sig_x = .1;
sig_e = sqrt(Model.getEvidenceVariance(like));
sig_e_C = sqrt(sig_x^2 + sig_e^2);

for si=1:length(sampleses)
    n_samples = sampleses(si);
    savename = sprintf('is_bias_%.2f_%.2f_%d_%d.mat', pri, like, n_samples, n_gridpts);
    savefile = fullfile(savedir, savename);
    if use_precomputed && exist(savefile, 'file')
        results = load(savefile);
        bias2_results{si} = results.bias.^2;
        variance_results{si} = results.variance;
    else
        bias = zeros(size(pp));
        variance = zeros(size(pp));
        parfor i=1:numel(pp)
            [true_update, sampled_updates] = ...
                run(pp(i), ee(i), gg(i), n_samples, n_repeats, sig_x, sig_e, pri);
            bias(i) = mean(sampled_updates) - true_update;
            variance(i) = var(sampled_updates, 1);
        end
        bias2_results{si} = bias.^2;
        variance_results{si} = variance;
        
        save(savefile, 'bias', 'variance');
    end
    
    %% Marginalize over pri/like distribution of e.
    
    prior_e = mog.create([+1, -1], [sig_e_C, sig_e_C], [pri, 1-pri]);
    pdf_e = mog.pdf(e_values, prior_e, true);
    
    % Note: meshgrid creates unintuitive slicing. Despite order of
    % arguments to meshgrid as [ee, pp, gg], ee varies over the second
    % dimension: assert(all(ee(1, :, 1) == e_values));
    marginal_bias2{si} = zeros(n_gridpts, n_gridpts);
    marginal_variance{si} = zeros(n_gridpts, n_gridpts);
    for i=1:n_gridpts
        for j=1:n_gridpts
            marginal_bias2{si}(i, j) = dot(squeeze(bias2_results{si}(i, :, j)), pdf_e);
            marginal_variance{si}(i, j) = dot(squeeze(variance_results{si}(i, :, j)), pdf_e);
        end
    end
end

%% Plot bias and variance as a function of gamma and num samples

figure;
for si=1:length(sampleses)
    subplot(1, 3, 1); hold on;
    plot(gamma_values, squeeze(mean(marginal_bias2{si})));
    title('Bias^2');
    
    subplot(1, 3, 2); hold on;
    plot(gamma_values, squeeze(mean(marginal_variance{si})));
    title('Variance');
    
    subplot(1, 3, 3); hold on;
    plot(gamma_values, squeeze(mean(marginal_bias2{si}) + mean(marginal_variance{si})));
    title('Bias^2 + Variance');
end

for sp=1:3
    subplot(1, 3, sp);
    xlabel('gamma');
    legend(arrayfun(@(s) [num2str(s) ' samples'], sampleses, 'UniformOutput', false), ...
        'Location', 'best');
end

end

function [true_update, sampled_updates] = run(lpo, e, gamma, n_samples, n_repeats, sig_x, sig_e, pri)

p_positive = exp(lpo) / (1 + exp(lpo));
prior_x = mog.create([-1, +1], [sig_x, sig_x], [1-p_positive p_positive]);
likelihood = mog.create(e, sig_e, 1);
Q = mog.prod(likelihood, prior_x);
sampled_updates = zeros(1, n_repeats);
for r=1:n_repeats
    xs = mog.sample(Q, n_samples);
    % Compute updates for C=+/-1
    updates_p = normpdf(xs, +1, sig_x);
    updates_m = normpdf(xs, -1, sig_x);
    % Compute importance-sampling weights: 1/prior.
    ws = 1 ./ mog.pdf(xs, prior_x, true);
    sampled_updates(r) = log(dot(ws, updates_p)) - log(dot(ws, updates_m));
end
sampled_updates = sampled_updates - gamma * lpo;
% true update based on p(e|C=+1), which is a mixture of gaussians with
% variance sig_x^2 + sig_e^2
sig_e_C = sqrt(sig_x^2 + sig_e^2);
p_e_C_p = mog.create([-1, +1], [sig_e_C, sig_e_C], [1-pri, pri]);
p_e_C_m = mog.create([-1, +1], [sig_e_C, sig_e_C], [pri, 1-pri]);
true_update = log(mog.pdf(e, p_e_C_p) / mog.pdf(e, p_e_C_m));
end