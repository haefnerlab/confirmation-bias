function [correct] = plotPriorLikelihoodSpace(priors, likelihoods, sampling_params)

% Get cartesian product of priors x likelihoods
[pp, ll] = meshgrid(priors, likelihoods);
correct = zeros(size(pp));
pp = pp(:);
ll = ll(:);

% Generate data
trials = 1000;
frames = 12;
prefix = sprintf('%dx%d_norminvlike', trials, frames);

parfor i=1:numel(pp)
    prior = pp(i);
    likelihood = ll(i);
    [data, var_e] = Model.genDataWithPriorLikelihood(trials, frames, prior, likelihood);
    sampling_params_copy = sampling_params;
    sampling_params_copy.var_e = var_e;
    sampling_params_copy.p_x = prior;
    results = Model.loadOrRunSamplingModel(data, prefix, sampling_params_copy);
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
saveas(gcf, 'PLSpace.fig');
end