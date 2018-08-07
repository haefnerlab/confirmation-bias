function plotTrueBiasCorrection(params, linespec)
savedir = fullfile('+Model', 'saved results');
if ~exist(savedir, 'dir'), mkdir(savedir); end

%% Get model bias for each combination of (LPO, evidence), then marginalize over evidence

lpo = linspace(-3, 3);
e_max = 1 + 5 * sqrt(params.var_s);
e_vals = linspace(-e_max, e_max, 1001)';

sig_s = sqrt(params.var_s);
data_pdf = mog.create([-1 +1], [sig_s sig_s], [params.prior_C 1-params.prior_C]);
e_pdf = mog.pdf(e_vals, data_pdf, true);

ideal_update = Model.logLikelihoodOdds(params, e_vals);
biases = zeros(length(lpo), length(e_vals));
for i=length(lpo):-1:1
    switch params.model
        case 'is'
            model_update = Model.isLogOddsUpdate(params, e_vals, lpo(i)*ones(size(e_vals)));
        case 'vb'
            model_update = Model.vbLogOddsUpdate(params, e_vals, lpo(i)*ones(size(e_vals)));
        case 'ideal'
            % Not sure why you would want to do this, but it's included for completeness
            model_update = ideal_update;
    end
    biases(i, :) = model_update - ideal_update;
end

expected_bias = biases * e_pdf;
% expected_bias = smooth(biases * e_pdf);
sigma_bias = sqrt((biases - expected_bias).^2 * e_pdf / sqrt(length(e_pdf)));

%% Get optimal value of linear gamma correction

uid = Model.getModelStringID(params);
optim_prefix = Model.getOptimPrefix({'gamma'}, 21);
optim_params = LoadOrRun(@Model.optimizeParams, {params, {'gamma'}, 21}, ...
    fullfile(savedir, [optim_prefix '_' uid]));

%% Plot it

if ~exist('linespec', 'var'), linespec = {'-'}; end

errorbar(lpo, expected_bias, sigma_bias, linespec{:});
xlabel('LPO');
ylabel('E_{model}[LLO] - E_{ideal}[LLO]');
title('True Nonlinear Bias Correction');

hold on;
h = plot(lpo, optim_params.gamma * lpo, 'k', 'LineStyle', '--', 'LineWidth', 2);
legend(h, {'best \gamma'});

end