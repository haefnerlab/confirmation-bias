function [correct, pk_ab, pk_tau] = plotPriorLikelihoodSpace(trials, frames, prior, likelihood, params, ideal_observer, pk)
%PLOTPRIORLIKELIHOODSPACE make a info_prior vs info_likelihood plot for the
%given params.
%
%params.p_match will take on values from 'prior'.
%params.var_e will take on values such that the AUC of x|e is 'likelihood'

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 6, ideal_observer = false; end
if nargin < 7, pk = false; end

[ll, pp] = meshgrid(likelihood, prior);

% Preallocate return values.
correct = nan(size(ll));
pk_ab = nan(size(ll));
pk_tau = nan(size(ll));

parfor i=1:numel(ll)
    params_copy = params;
    % Set variances for this pair of prior & likelihood values.
    params_copy.var_e = Model.getEvidenceVariance(ll(i));
    params_copy.p_match = pp(i);
    
    % Generate data and run the model.
    [data, prefix] = Model.genDataWithParams(trials, frames, params_copy);
    if ~ideal_observer
        results = Model.loadOrRunSamplingModel(data, prefix, params_copy);
    else
        results = Model.runIdealObserver(data, params_copy);
    end
    % Fit PK if requested
    if pk
        [~, ~, expfit] = Model.loadOrRunModelPK(Model.getModelStringID(prefix, params_copy), data, results, [1 0 10]);
        pk_ab(i) = expfit(1) / expfit(2);
        pk_tau(i) = 1 / expfit(3);
    end
    correct(i) = sum(results.choices == +1) / trials;
end

% Plot percent correct
figure();
imagesc(correct, [0.5 1.0]); axis image; colorbar;
% Add contour line at threshold
hold on; contour(imfilter(correct, smoothkernel(5), 'replicate'), [0.7 0.7], '-w', 'LineWidth', 2);
prior_tick_indices = round(linspace(1, length(prior), min(length(prior), 5)));
like_tick_indices = round(linspace(1, length(likelihood), min(length(likelihood), 5)));
set(gca, 'YTick', prior_tick_indices);
set(gca, 'XTick', like_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, prior(prior_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, likelihood(like_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('P_{likelihood}');
ylabel('P_{prior}');
title('PL-Space');
if ~ideal_observer
    figname = sprintf('PLSpace_%dx%d_vx%.2f_pD%.2f_gam%.2f_ns%d_b%d.fig', ...
        trials, frames, params.var_x, params.prior_D, ...
        params.gamma, params.samples, params.batch);
else
    figname = sprintf('PLSpace_%dx%d_vx%.2f_ideal.fig', trials, frames, params.var_x);
end

saveas(gcf, fullfile(savedir, figname));

% Plot PK fit terms.
if pk
    figure();
    subplot(1,2,1);
    imagesc(pk_ab, [0 10]); axis image; colorbar;
    set(gca, 'YTick', prior_tick_indices);
    set(gca, 'XTick', like_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, prior(prior_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, likelihood(like_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('P_{likelihood}');
    ylabel('P_{prior}');
    title('PL-Space: A/B of PK fit');
    
    subplot(1,2,2);
    imagesc(pk_tau, [0 5]); axis image; colorbar;
    set(gca, 'YTick', prior_tick_indices);
    set(gca, 'XTick', like_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, prior(prior_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, likelihood(like_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('P_{likelihood}');
    ylabel('P_{prior}');
    title('PL-Space: 1/\tau of PK fit');
    
    if ~ideal_observer
        figname = sprintf('PLSpace_PK_%dx%d_vx%.2f_pD%.2f_gam%.2f_ns%d_b%d.fig', ...
            trials, frames, params.var_x, params.prior_D, ...
            params.gamma, params.samples, params.batch);
    else
        figname = sprintf('PLSpace_PK_%dx%d_vx%.2f_ideal.fig', trials, frames, params.var_x);
    end
    saveas(gcf, fullfile(savedir, figname));
end
end

function kernel = smoothkernel(n)
x = linspace(-2,2,n);
[xx,yy] = meshgrid(x,x);
kernel = exp(-(xx.^2 + yy.^2));
kernel = kernel / sum(kernel(:));
end