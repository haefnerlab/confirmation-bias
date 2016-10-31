function [correct] = plotPriorLikelihoodSpace(priors, likelihoods, sampling_params)

% Get cartesian product of priors x likelihoods
[pp, ll] = meshgrid(priors, likelihoods);
correct = zeros(size(pp));
pp = pp(:);
ll = ll(:);

% Generate data
trials = 1000;
frames = 12;
prefix = sprintf('%dx%d_hackedlike', trials, frames);

parfor i=1:numel(pp)
    prior = pp(i);
    likelihood = ll(i);
    [data, var_e] = gen_data(trials, frames, prior, likelihood);
    sampling_params_copy = sampling_params;
    sampling_params_copy.var_e = var_e;
    sampling_params_copy.p_x = prior;
    results = Model.loadOrRunSamplingModel(data, prefix, sampling_params_copy);
    correct(i) = sum(results.choices == +1) / trials;
end

imagesc(correct); axis image; colorbar;
set(gca, 'XTickLabel', arrayfun(@num2str, likelihoods, 'UniformOutput', false));
set(gca, 'YTickLabel', arrayfun(@num2str, priors, 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('P_{likelihood}');
ylabel('P_{prior}');
title('Percent Correct in prior-likelihood space');
saveas(gcf, 'PLSpace.fig');
end

function [data, var_e] = gen_data(trials, frames, prior, likelihood)
%GEN_DATA generates a set of trials, all with correct choice +1
data = zeros(trials, frames);
var_e = 1 / (likelihood - 0.5); % TODO
stdev_e = sqrt(var_e);
for t=1:trials
    % generate the 'center' of each frame according to 'prior'
    centers = stdev_e * sign(prior - rand(1, frames));
    % draw signal from around the center with variance set by 'likelihood'
    data(t, :) = normrnd(centers, stdev_e);
end
end