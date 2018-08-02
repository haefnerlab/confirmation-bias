function plotCSSlopes(category_infos, sensory_infos, params, beta_range)

savedir = fullfile('+VariationalModel', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

[ss, cc] = meshgrid(sensory_infos, category_infos);

% Preallocate return values.
corrects = nan(size(ss));
slopes = nan(size(ss));
slopeErrors = nan(size(ss));

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
    
    % Run the model
    results_uid = VariationalModel.getModelStringID(params_copy);
    results = LoadOrRun(params.model_fun, {params_copy}, fullfile(params.save_dir, results_uid));
    
    data = SamplingModel.genDataWithParams(results.params);
    [data, choices] = flipTrials(data, results.choices);
    weights = CustomRegression.PsychophysicalKernel(data, choices, 0, 0, 0);
    [expfit, expErrors] = CustomRegression.expFit(weights);
    corrects(i) = mean(results.choices == +1);
    slopes(i) = expfit(2);
    slopeErrors(i) = expErrors(2);
end

colors = arrayfun(@(b) SamplingModel.betacolor(b, beta_range(1), beta_range(2)), linspace(beta_range(1), beta_range(2)), ...
    'UniformOutput', false);
cmap = vertcat(colors{:});

% gray_alphas = 1 - exp(-slopeErrors.^2/6);
gray_alphas = 1 - (1 + exp(-(corrects-.55)*100)).^-1;
gray_lum = 1;
gray_rgbs = gray_lum * ones([size(ss) 3]);

% Plot slope
fig = figure(); hold on;
% Slopes contourf
[c,h] = contourf(smoothn(slopes), linspace(beta_range(1), beta_range(2), 11));
axis image;
colormap(cmap);
clabel(c,h,'Color',[1 1 1]);
% Gray overlay
h = image(gray_rgbs);
set(h, 'AlphaData', gray_alphas);
% Labels n such
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
fullname = [VariationalModel.getModelStringID(params) '.fig'];
figname = regexprep(fullname, '_cinfo[\d.]+_sinfo[\d.]+', '');
saveas(fig, fullfile(savedir, figname));
end

function [data, choices] = flipTrials(data, choices)
flip_indexes = rand(length(choices), 1) < 0.5;
data(flip_indexes, :) = -data(flip_indexes, :);
choices = choices == +1;
choices(flip_indexes) = ~choices(flip_indexes);
end