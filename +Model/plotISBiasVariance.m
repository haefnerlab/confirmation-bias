function [marg_bias2_results, marg_variance_results, e_values, log_prior_c, gamma_values] = ...
    plotISBiasVariance(cat_info, sense_info, sampleses, n_gridpts, n_repeats, use_precomputed, verbose)

if nargin < 6, use_precomputed = false; end
if nargin < 7, verbose = false; end

%% Load or compute results
savedir = fullfile('+Model', 'saved results');

log_prior_c = linspace(-2, 2, n_gridpts);
e_values = linspace(-2, 2, n_gridpts);
gamma_values = linspace(0, 1, n_gridpts);
[ee, pp, gg] = meshgrid(e_values, log_prior_c, gamma_values);

% Results will be cat points x sense points x samples x gammas
res_size = [length(cat_info) length(sense_info) length(sampleses) n_gridpts];
marg_bias2_results = zeros(res_size);
marg_variance_results = zeros(res_size);

sig_x = .1;

for cat_idx=1:length(cat_info)
    cat = cat_info(cat_idx);
    for sense_idx=1:length(sense_info)
        sens = sense_info(sense_idx);
        sig_e = sqrt(Model.getEvidenceVariance(sens));
        sig_e_C = sqrt(sig_x^2 + sig_e^2);
        for samp_idx=1:length(sampleses)
            n_samples = sampleses(samp_idx);
            savename = sprintf('is_bias_%.03f_%.03f_%d_%d.mat', ...
                cat, sens, n_samples, n_gridpts);
            savefile = fullfile(savedir, savename);
            if use_precomputed && exist(savefile, 'file')
                if verbose, disp(['Loading ' savename]); end
                load(savefile);
            else
                if verbose, disp(['Computing ' savename]); end
                bias = zeros(size(pp(:,:,1)));
                variance = zeros(size(bias));
                parfor i=1:numel(ee(:,:,1))
                    % In parallel get updates for e/lpo combinations. Gamma
                    % term added later in vectorized form.
                    [true_update, sampled_updates] = ...
                        run(pp(i), ee(i), n_samples, n_repeats, sig_x, sig_e, cat);
                    bias(i) = mean(sampled_updates) - true_update;
                    variance(i) = var(sampled_updates, 1);
                end
                bias = repmat(bias, [1 1 length(gamma_values)]) - pp .* gg;
                variance = repmat(variance, [1 1 length(gamma_values)]);
                save(savefile, 'bias', 'variance');
            end
            [marg_bias2_results(cat_idx, sense_idx, samp_idx, :), marg_variance_results(cat_idx, sense_idx, samp_idx, :)] = ...
                get_marginal_results(bias.^2, variance, sig_e_C);
        end
    end
end

%% Plot MSE for low/med/high gamma x low/med/high num samples
bfig = figure;
vfig = figure;
mfig = figure;

idx_samples = round([1, length(sampleses)/2, length(sampleses)]);
idx_gammas = round([1, n_gridpts/2, n_gridpts]);

marg_mse = marg_bias2_results + marg_variance_results;
bias2_range = [0 max(marg_bias2_results(:))];
variance_range = [0 max(marg_variance_results(:))];
mse_range = [0 max(marg_mse(:))];

for i=1:3
    for j=1:3
        figure(bfig);
        subplot(3, 3, sub2ind([3 3], i, j));
        imagesc(marg_bias2_results(:, :, i, j), bias2_range);
        colorbar;
        axis image; set(gca, 'YDir', 'normal');
        title(sprintf('Bias^2 (%d samples, %.2f gamma)', sampleses(idx_samples(i)), gamma_values(idx_gammas(j))));
        xlabel('Sensory Info');
        set(gca, 'XTick', 1:2:length(sense_info));
        set(gca, 'XTickLabel', sense_info(1:2:length(sense_info)));
        ylabel('Category Info');
        set(gca, 'YTick', 1:2:length(cat_info));
        set(gca, 'YTickLabel', cat_info(1:2:length(cat_info)));

        figure(vfig);
        subplot(3, 3, sub2ind([3 3], i, j));
        imagesc(marg_variance_results(:, :, i, j), variance_range);
        colorbar;
        axis image; set(gca, 'YDir', 'normal');
        title(sprintf('Variance (%d samples, %.2f gamma)', sampleses(idx_samples(i)), gamma_values(idx_gammas(j))));
        xlabel('Sensory Info');
        set(gca, 'XTick', 1:2:length(sense_info));
        set(gca, 'XTickLabel', sense_info(1:2:length(sense_info)));
        ylabel('Category Info');
        set(gca, 'YTick', 1:2:length(cat_info));
        set(gca, 'YTickLabel', cat_info(1:2:length(cat_info)));

        figure(mfig);
        subplot(3, 3, sub2ind([3 3], i, j));
        imagesc(marg_mse(:, :, i, j), mse_range);
        colorbar;
        axis image; set(gca, 'YDir', 'normal');
        title(sprintf('MSE (%d samples, %.2f gamma)', sampleses(idx_samples(i)), gamma_values(idx_gammas(j))));
        xlabel('Sensory Info');
        set(gca, 'XTick', 1:2:length(sense_info));
        set(gca, 'XTickLabel', sense_info(1:2:length(sense_info)));
        ylabel('Category Info');
        set(gca, 'YTick', 1:2:length(cat_info));
        set(gca, 'YTickLabel', cat_info(1:2:length(cat_info)));
    end
end

end

function [marg_bias2, marg_variance] = get_marginal_results(bias2, variance, sig_e_C)
% Marginalize results over cat/sense distribution of e, and take mean over
% possible values of LPO. Returns one result per 'gamma' value.

n_gridpts = size(bias2, 1);
marg_bias2 = zeros(1, n_gridpts);
marg_variance = zeros(1, n_gridpts);

e_values = linspace(-2, 2, n_gridpts);
prior_e = mog.create([+1, -1], [sig_e_C, sig_e_C], [.5 .5]);
pdf_e = mog.pdf(e_values, prior_e, true);

% Note: meshgrid creates unintuitive slicing. Despite order of
% arguments to meshgrid as [ee, pp, gg], ee varies over the second
% dimension: above, we could do `assert(all(ee(1, :, 1) == e_values));`
for gam_idx=1:n_gridpts
    bb = zeros(1, n_gridpts);
    vv = zeros(1, n_gridpts);
    for lpo_idx=1:n_gridpts
        bb(lpo_idx) = dot(squeeze(bias2(lpo_idx, :, gam_idx)), pdf_e);
        vv(lpo_idx) = dot(squeeze(variance(lpo_idx, :, gam_idx)), pdf_e);
    end
    % Simple mean over LPO dimension
    % TODO - what is the distribution of LPO values here?
    marg_bias2(gam_idx) = mean(bb);
    marg_variance(gam_idx) = mean(vv);
end

end

function [true_update, sampled_updates] = run(lpo, e, n_samples, n_repeats, sig_x, sig_e, cat_info)
p_positive = exp(lpo) / (1 + exp(lpo));
prior_x = mog.create([-1, +1], [sig_x, sig_x], [1-p_positive p_positive]);
likelihood = mog.create(e, sig_e, 1);
Q = mog.prod(likelihood, prior_x);
sampled_updates = zeros(1, n_repeats);
p_x_Cp = mog.create([-1, +1], [sig_x, sig_x], [1-cat_info, cat_info]);
p_x_Cm = mog.create([-1, +1], [sig_x, sig_x], [cat_info, 1-cat_info]);
for r=1:n_repeats
    xs = mog.sample(Q, n_samples);
    % Compute updates for C=+/-1 based on p(x|C+) and p(x|C-)
    updates_p = mog.pdf(xs, p_x_Cp);
    updates_m = mog.pdf(xs, p_x_Cm);
    % Compute importance-sampling weights: 1/prior.
    ws = 1 ./ mog.pdf(xs, prior_x);
    sampled_updates(r) = log(dot(ws, updates_p)) - log(dot(ws, updates_m));
end
% true update based on p(e|C=+1), which is a mixture of gaussians with
% variance sig_x^2 + sig_e^2
sig_e_C = sqrt(sig_x^2 + sig_e^2);
p_e_Cp = mog.create([-1, +1], [sig_e_C, sig_e_C], [1-cat_info, cat_info]);
p_e_Cm = mog.create([-1, +1], [sig_e_C, sig_e_C], [cat_info, 1-cat_info]);
true_update = mog.logpdf(e, p_e_Cp) - mog.logpdf(e, p_e_Cm);
end