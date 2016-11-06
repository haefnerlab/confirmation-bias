function [correct] = plotPriorLikelihoodSpace(trials, frames, I_prior, I_likelihood, params, ideal_observer, pk)
%PLOTPRIORLIKELIHOODSPACE make a info_prior vs info_likelihood plot for the
%given params.
%
%params.p_match and params.var_e will be set based on alpha and beta.

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 6, ideal_observer = false; end
if nargin < 7, pk = false; end

% Check for bounds on alpha and beta imposed by var_x
[alpha_min, beta_max] = getInformationBounds(params.var_x);

% Adjust I_prior range to fit; assume logspace.
alpha = max(min(I_prior), alpha_min);
beta = min(max(I_prior), beta_max);
I_prior = logspace(log(alpha)/log(10), log(beta)/log(10), length(I_prior));

% Get meshgrid of alpha:beta for both the prior (y) and likelihood (x)
% axes.
[I_L, I_P] = meshgrid(I_likelihood, I_prior);

correct = nan(size(I_P));
pk_ab = zeros(size(I_P));
pk_tau = zeros(size(I_P));

parfor i=1:numel(I_P)
    params_copy = params;
    % Set variances for given 'information' contents.
    params_copy.var_e = 1 / I_L(i);
    params_copy.p_match = solveForPriorParams(I_P(i), params.var_x);
    
    % Generate data and run the model.
    [data, prefix] = Model.genDataWithParams(trials, frames, params_copy);
    if ~ideal_observer
        results = Model.loadOrRunSamplingModel(data, prefix, params_copy);
    else
        results = Model.runIdealObserver(data, params_copy);
    end
    % Fit PK if requested
    if pk
        [weights, ~] = Model.loadOrRunModelPK(Model.getModelStringID(prefix, params_copy), data, results, [0 0 10]);
        expfit = CustomRegression.expFit(weights);
        pk_ab(i) = expfit(1) / expfit(2);
        pk_tau(i) = 1 / expfit(3);
    end
    correct(i) = sum(results.choices == +1) / trials;
end

% Plot percent correct
figure();
imagesc(correct, [0.5 1.0]); axis image; colorbar;
prior_tick_indices = round(linspace(1, length(I_prior), min(length(I_prior), 5)));
like_tick_indices = round(linspace(1, length(I_likelihood), min(length(I_likelihood), 5)));
set(gca, 'YTick', prior_tick_indices);
set(gca, 'XTick', like_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, I_prior(prior_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, I_likelihood(like_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('Info_{likelihood}');
ylabel('Info_{prior}');
title('PL-Space');
if ~ideal_observer
    figname = sprintf('PLSpace_%dx%d_vx%.2f_pD%.2f_gam%.2f_S%d.fig', ...
        trials, frames, params.var_x, params.prior_D, ...
        params.gamma, params.samples);
else
    figname = sprintf('PLSpace_%dx%d_vx%.2f_ideal.fig', trials, frames, params.var_x);
end

saveas(gcf, fullfile(savedir, figname));

% Plot PK fit terms.
if pk
    figure();
    subplot(1,2,1);
    imagesc(pk_ab); axis image; colorbar;
    set(gca, 'YTick', prior_tick_indices);
    set(gca, 'XTick', like_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, I_prior(prior_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, I_likelihood(like_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('Info_{likelihood}');
    ylabel('Info_{prior}');
    title('PL-Space: A/B of PK fit');

    subplot(1,2,2);
    imagesc(pk_tau); axis image; colorbar;
    set(gca, 'YTick', prior_tick_indices);
    set(gca, 'XTick', like_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, I_prior(prior_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, I_likelihood(like_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('Info_{likelihood}');
    ylabel('Info_{prior}');
    title('PL-Space: 1/\tau of PK fit');

    if ~ideal_observer
        figname = sprintf('PLSpace_PK_%dx%d_vx%.2f_pD%.2f_gam%.2f_S%d.fig', ...
            trials, frames, params.var_x, params.prior_D, ...
            params.gamma, params.samples);
    else
        figname = sprintf('PLSpace_PK_%dx%d_vx%.2f_ideal.fig', trials, frames, params.var_x);
    end
    saveas(gcf, fullfile(savedir, figname));
end
end

function [alpha_min, beta_max] = getInformationBounds(var_x)
% There is a limit to the variance that can be contributed by p_match at a
% given value of var_x. This function returns the low and high bounds on
% information.
alpha_min = 1 / (var_x + 0.5); % maximal variance of 2*0.25 when p is 0.5.
beta_max = 1 / var_x; % minimal variance of var_x when p is 1.
end

function p_match = solveForPriorParams(info_prior, var_x)
p_match = real((1 + sqrt(1 + 2*(var_x - 1./info_prior))) / 2);
p_match = min(max(p_match, 0.5), 1);
end