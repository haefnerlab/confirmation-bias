function [slopes, slopeErrors, corrects] = plotCSSlopes(category_infos, sensory_infos, params, optimize, optim_grid_size, rb_range)

savedir = fullfile('+SamplingModel', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end
if nargin < 4, optimize = {}; end
if nargin < 5, optim_grid_size = 11; end

[ss, cc] = meshgrid(sensory_infos, category_infos);

% Preallocate return values.
corrects = nan(size(ss));
slopes = nan(size(ss));
slopeErrors = nan(size(ss));

optim_prefix = SamplingModel.getOptimPrefix(optimize, optim_grid_size);

parfor i=1:numel(ss)
    params_copy = params;
    % Set data-generating parameters.
    params_copy.sensory_info = ss(i);
    params_copy.category_info = cc(i);
    % Set variances for this pair of category- & sensory-info values. (That is, assume that the
    % model 'knows' the environment statistics)
    params_copy.var_s = SamplingModel.getEvidenceVariance(ss(i));
    params_copy.p_match = cc(i);
    
    % TODO - smarter setting of seed?
    params_copy.seed = randi(1000000000);
    
    % Run the model if needed
    results_uid = SamplingModel.getModelStringID(params_copy);
    
    if isempty(optimize)
        results = LoadOrRun(@SamplingModel.runSamplingModelFast, {params_copy}, ...
            fullfile(params.save_dir, results_uid));
    else
        % Find optimal param settings.
        [optim_params, ~] = LoadOrRun(@SamplingModel.optimizeParams, ...
            {params_copy, optimize, optim_grid_size}, ...
            fullfile(params.save_dir, [optim_prefix '_' results_uid]));
        % Get model results at the optimal param settings.
        best_results_uid = SamplingModel.getModelStringID(optim_params);
        results = LoadOrRun(@SamplingModel.runSamplingModelFast, {optim_params}, ...
            fullfile(params.save_dir, best_results_uid));
    end
    
    data = SamplingModel.genDataWithParams(results.params);
    [data, choices] = flipTrials(data, results.choices);
    weights = CustomRegression.PsychophysicalKernel(data, choices, 0, 0, 0);
    [expfit, expErrors] = CustomRegression.expFit(weights);
    corrects(i) = mean(results.choices == +1);
    slopes(i) = expfit(2);
    slopeErrors(i) = expErrors(2);
end

colors = arrayfun(@(b) SamplingModel.betacolor(b, rb_range(1), rb_range(2)), linspace(rb_range(1), rb_range(2)), ...
    'UniformOutput', false);
cmap = vertcat(colors{:});

% gray_alphas = 1 - exp(-slopeErrors.^2/6);
gray_alphas = 1 - (1 + exp(-(corrects-.55)*100)).^-1;
gray_lum = 1;
gray_rgbs = gray_lum * ones([size(ss) 3]);

% Plot slope
fig = figure(); hold on;
[c,h] = contourf(smoothn(slopes), linspace(rb_range(1), rb_range(2), 11));
axis image;
colormap(cmap);
clabel(c,h);
% imagesc(slopes, range); axis image; colormap(cmap); colorbar;
% h = image(gray_rgbs);
% set(h, 'AlphaData', gray_alphas);
category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
set(gca, 'YTick', category_tick_indices);
set(gca, 'XTick', sensory_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('Sensory Info');
ylabel('Category Info');
title('C-S Space: PK Slope (\beta)');
figname = sprintf('CSSpace_slope_%dx%d_%s_vx%.2f_pC%.2f_gam%.2f_ns%d_b%d_%d_%.2e.fig', ...
    params.trials, params.frames, optim_prefix, params.var_x, params.prior_C, ...
    params.gamma, params.samples, params.batch, params.importance_norm, params.noise);
saveas(fig, fullfile(savedir, figname));
end

function [data, choices] = flipTrials(data, choices)
flip_indexes = rand(length(choices), 1) < 0.5;
data(flip_indexes, :) = -data(flip_indexes, :);
choices = choices == +1;
choices(flip_indexes) = ~choices(flip_indexes);
end