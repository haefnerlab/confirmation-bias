function [correct, fig, opt_fig] = plotCategorySensorySpace(category_infos, sensory_infos, params, optimize, optim_grid_size, do_plot)
%PLOTCATEGORYSENSORYSPACE make category_info vs sensory_info plots for the
%given params.

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 4, optimize = {}; end
if nargin < 5, optim_grid_size = 11; end
if nargin < 6, do_plot = true; end

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
    params_copy = Model.setCategorySensoryInfo(params, cc(i), ss(i));
    
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

if do_plot
    % Plot percent correct
    fig = figure();
    imagesc('XData', sensory_infos, 'YData', category_infos, 'CData', correct, [0.5 1.0]);
    axis equal; axis image; colorbar;
    % Add contour line at threshold
    [~, threshold_pts] = Model.getThresholdPoints(category_infos, params, .7, 50);
    hold on; plot(threshold_pts(:,1), threshold_pts(:,2), '-w', 'LineWidth', 2);
    category_tick = category_infos(round(linspace(1, length(category_infos), 5)));
    sensory_tick = sensory_infos(round(linspace(1, length(sensory_infos), 5)));
    set(gca, 'YTick', category_tick);
    set(gca, 'XTick', sensory_tick);
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
        category_tick = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
        sensory_tick = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
        set(gca, 'YTick', category_tick);
        set(gca, 'XTick', sensory_tick);
        set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick), 'UniformOutput', false));
        set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick), 'UniformOutput', false));
        set(gca, 'YDir', 'Normal');
        xlabel('SI');
        ylabel('CI');
        title(['Optimized value of ' optimize{i}]);
        figname = ['CSSpace_optim_' optimize{i} '_' Model.getModelStringID(params, true) '.fig'];
        saveas(gcf, fullfile(savedir, figname));
    end
end
end