function [correct, fig, opt_fig] = plotCategorySensorySpace(category_infos, sensory_infos, params, optimize, optim_grid_size)
%PLOTCATEGORYSENSORYSPACE make category_info vs sensory_info plots for the
%given params.

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 4, optimize = {}; end
if nargin < 5, optim_grid_size = 11; end

if exist('optimize', 'var') && ~isempty(optimize)
    optim_prefix = Model.getOptimPrefix(optimize, optim_grid_size);
    results_uid = Model.getModelStringID(params, true);
    optim_results_uid = ['optim_CS_' optim_prefix '_' results_uid];
    [optim_params, ~, ~] = LoadOrRun(@Model.optimizeParamsOverCS, ...
        {category_infos, sensory_infos, params, optimize, optim_grid_size}, ...
        fullfile(params.save_dir, optim_results_uid));
else
    optimize = {};
    optim_params = [];
end

if strcmpi(params.model, 'ideal') && ~isempty(optimize)
    error('Nothing to optimize for the ideal observer');
end

[ss, cc] = meshgrid(sensory_infos, category_infos);

% Preallocate return values.
correct = nan(size(ss));

for i=1:numel(ss)
    params_copy = params;
    % Set data-generating parameters.
    params_copy.sensory_info = ss(i);
    params_copy.category_info = cc(i);
    % Set variances for this pair of category- & sensory-info values. (That is, assume that the
    % model 'knows' the environment statistics)
    params_copy.var_s = Model.getEvidenceVariance(ss(i));
    params_copy.p_match = cc(i);
    
    % TODO - smarter setting of seed?
    params_copy.seed = randi(1000000000);
    
    if isempty(optimize)
        % Run the model on the given params.
        results_uid = Model.getModelStringID(params_copy);
        results = LoadOrRun(@Model.runVectorized, {params_copy}, ...
            fullfile(params.save_dir, results_uid));
    else
        % Run the model using the best params at this value of category and sensory info.
        for iVar=1:length(optimize)
            params_copy.(optimize{iVar}) = optim_params(i).(optimize{iVar});
        end
        results_uid = Model.getModelStringID(params_copy);
        results = LoadOrRun(@Model.runVectorized, {params_copy}, ...
            fullfile(params.save_dir, results_uid), '-verbose');
    end
    
    [~, correct_categories] = Model.genDataWithParams(results.params);
    
    correct(i) = mean(results.choices == correct_categories);
    if mod(i,10)==1, disp(i); end
end

% Plot percent correct
fig = figure();
imagesc(correct, [0.5 1.0]); axis image; colorbar;
% Add contour line at threshold
hold on; contour(imfilter(correct, smoothkernel(9), 'replicate'), [0.7 0.7], '-w', 'LineWidth', 2);
category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
set(gca, 'YTick', category_tick_indices);
set(gca, 'XTick', sensory_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('SI');
ylabel('CI');
title('Percent Correct');

figname = ['CSSpace_' Model.getModelStringID(params, true) '.fig'];

saveas(gcf, fullfile(savedir, figname));

% Plot value of optimized parameters.
for i=length(optimize):-1:1
    opt_fig(i) = figure();
    
    % Unravel optimal param values
    opt_param_values = reshape([optim_params.(optimize{i})]', size(ss));
    
    imagesc(opt_param_values); axis image; colorbar; colormap cool;
    category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
    sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
    set(gca, 'YTick', category_tick_indices);
    set(gca, 'XTick', sensory_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('SI');
    ylabel('CI');
    title(['Optimized value of ' optimize{i}]);
    figname = ['CSSpace_optim_' optimize{i} '_' Model.getModelStringID(params, true) '.fig'];
    saveas(gcf, fullfile(savedir, figname));
end
end

function kernel = smoothkernel(n)
x = linspace(-2,2,n);
[xx,yy] = meshgrid(x,x);
kernel = exp(-(xx.^2 + yy.^2));
kernel = kernel / sum(kernel(:));
end