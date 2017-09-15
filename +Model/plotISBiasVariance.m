function plotISBiasVariance(n_samples, n_gridpts, n_repeats, use_precomputed)

if nargin < 4, use_precomputed = false; end

%% Load or compute results

log_prior_c = linspace(-2, 2, n_gridpts);
e_values = linspace(-.5, .5, n_gridpts);
gamma_values = linspace(0, 1, n_gridpts);

savedir = fullfile('+Model', 'saved results');
savename = sprintf('is_bias_%d_%d.mat', n_samples, n_gridpts);
savefile = fullfile(savedir, savename);
if use_precomputed && exist(savefile, 'file')
    results = load(savefile);
    bias_results = results.bias_results;
    variance_results = results.variance_results;
else
    
    [ee, pp, gg] = meshgrid(e_values, log_prior_c, gamma_values);
    bias_results = zeros(size(pp));
    variance_results = zeros(size(pp));
    parfor i=1:numel(pp)
        [true_update, sampled_updates] = ...
            run(pp(i), ee(i), gg(i), n_samples, n_repeats);
        bias_results(i) = mean(sampled_updates) - true_update;
        variance_results(i) = var(sampled_updates, 1);
    end
    
    save(savefile, 'bias_results', 'variance_results');
end

%% Marginalize out e according to some experimenter distribution around 0.

gen_sig_e = .1;
pdf_e = normpdf(e_values, 0, gen_sig_e);
pdf_e = pdf_e / sum(pdf_e);

marginal_bias = sum(bias_results .* pdf_e(:), 1);
marginal_variance = sum(variance_results .* pdf_e(:), 1);

%% Plot bias and variance images
tick_idxs = round([1, n_gridpts/4, n_gridpts/2, n_gridpts/4*3, n_gridpts]);

figure;
subplot(1, 2, 1);
imagesc(marginal_bias);
imagesc(marginal_variance);
xlabel('log prior odds');
ylabel('gamma');
title([num2str(n_samples) ' samples']);

set(gca, 'XTick', tick_idxs);
set(gca, 'XTickLabel', log_prior_c(tick_idxs));
set(gca, 'YTick', tick_idxs);
set(gca, 'YTickLabel', e_values(tick_idxs));

%% Plot bias and variance vs log_prior_c

figure;
subplot(2, 2, 1);
imagesc(bias_results(:, :, 1).^2);
axis image;
colorbar;
xlabel('log P_{t-1}(C+)/P_{t-1}(C-)');
ylabel('e value');
title(['Bias^2 (' num2str(n_samples) ' samples, \gamma=0)']);

set(gca, 'XTick', tick_idxs);
set(gca, 'XTickLabel', log_prior_c(tick_idxs));
set(gca, 'YTick', tick_idxs);
set(gca, 'YTickLabel', e_values(tick_idxs));

subplot(2, 2, 3);
imagesc(bias_results(:, :, end).^2);
axis image;
colorbar;
xlabel('log P_{t-1}(C+)/P_{t-1}(C-)');
ylabel('e value');
title(['Bias^2 (' num2str(n_samples) ' samples, \gamma=1)']);

set(gca, 'XTick', tick_idxs);
set(gca, 'XTickLabel', log_prior_c(tick_idxs));
set(gca, 'YTick', tick_idxs);
set(gca, 'YTickLabel', e_values(tick_idxs));

subplot(2, 2, 2);
imagesc(variance_results(:, :, 1));
axis image;
colorbar;
xlabel('log P_{t-1}(C+)/P_{t-1}(C-)');
ylabel('e value');
title(['Variance (' num2str(n_samples) ' samples, \gamma=0)']);

set(gca, 'XTick', tick_idxs);
set(gca, 'XTickLabel', log_prior_c(tick_idxs));
set(gca, 'YTick', tick_idxs);
set(gca, 'YTickLabel', e_values(tick_idxs));

end

function [true_update, sampled_updates] = run(lpo, e, gamma, n_samples, n_repeats)
sig_x = .1;
sig_e = .5;

p_positive = exp(lpo) / (1 + exp(lpo));
prior = mog.create([-1, +1], [sig_x, sig_x], [1-p_positive p_positive]);
likelihood = mog.create(e, sig_e, 1);
Q = mog.prod(likelihood, prior);
sampled_updates = zeros(1, n_repeats);
for r=1:n_repeats
    xs = mog.sample(Q, n_samples);
    % Compute updates for C=+/-1
    updates_p = normpdf(xs, +1, sig_x);
    updates_m = normpdf(xs, -1, sig_x);
    % Compute importance-sampling weights: 1/prior.
    ws = 1 ./ mog.pdf(xs, prior, true);
    sampled_updates(r) = log(dot(ws, updates_p)) - log(dot(ws, updates_m));
end
% true update based on p(e|C=+1), which is a normal with variance
% sig_x^2 + sig_e^2
log_p_e_C_p = -0.5 * (+1 - e)^2 / (sig_x^2 + sig_e^2);
log_p_e_C_m = -0.5 * (-1 - e)^2 / (sig_x^2 + sig_e^2);
true_update = log_p_e_C_p - log_p_e_C_m;
sampled_updates = sampled_updates - gamma * lpo;
end