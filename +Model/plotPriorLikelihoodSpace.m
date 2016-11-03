function [correct] = plotPriorLikelihoodSpace(trials, frames, priors, likelihoods, sampling_params, ideal_observer)

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 6, ideal_observer = false; end

% Get cartesian product of priors x likelihoods
[ll, pp] = meshgrid(likelihoods, priors);
correct = zeros(size(pp));
pp = pp(:);
ll = ll(:);

parfor i=1:numel(pp)
    prior = pp(i);
    likelihood = ll(i);
    [data, var_e, prefix] = Model.genDataWithPriorLikelihood(trials, frames, prior, likelihood);
    sampling_params_copy = sampling_params;
    sampling_params_copy.var_e = var_e;
    sampling_params_copy.p_match = prior;
    if ~ideal_observer
        results = Model.loadOrRunSamplingModel(data, prefix, sampling_params_copy);
    else
        results = Model.runIdealObserver(data, sampling_params_copy);
    end
    correct(i) = sum(results.choices == +1) / trials;
end

sub_prior = 1:5:length(priors);
sub_like = 1:5:length(likelihoods);

imagesc(correct, [0.5 1.0]); axis image; colorbar;
set(gca, 'XTick', sub_like);
set(gca, 'YTick', sub_prior);
set(gca, 'XTickLabel', arrayfun(@num2str, likelihoods(sub_like), 'UniformOutput', false));
set(gca, 'YTickLabel', arrayfun(@num2str, priors(sub_prior), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('P_{likelihood}');
ylabel('P_{prior}');
title('Percent Correct in prior-likelihood space');
if ~ideal_observer
    figname = sprintf('PLSpace_%dx%d_vx%.2f_pD%.2f_gam%.2f_S%d.fig', ...
        trials, frames, sampling_params.var_x, sampling_params.prior_D, ...
        sampling_params.gamma, sampling_params.samples);
else
    figname = sprintf('PLSpace_%dx%d_ideal.fig', trials, frames);
end
saveas(gcf, fullfile(savedir, figname));
end