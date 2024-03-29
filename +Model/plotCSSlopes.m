function [slopes, slopeErrors, smoothSlopes, corrects, fig, cmap] = plotCSSlopes(category_infos, sensory_infos, params, rb_range, contour_pc, contour_through_sc, optim, optim_grid_size)

memodir = fullfile('+Model', 'saved results');
savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if exist('optim', 'var') && ~isempty(optim)
    optim_prefix = Model.getOptimPrefix(optim, optim_grid_size);
    results_uid = Model.getModelStringID(params, true);
    optim_results_uid = ['optim_CS_' optim_prefix '_' results_uid];
    [optim_params, ~] = LoadOrRun(@Model.optimizeParamsOverCS, ...
        {category_infos, sensory_infos, params, optim, optim_grid_size}, ...
        fullfile(params.save_dir, optim_results_uid));
else
    optim = {};
    optim_params = [];
    optim_prefix = '';
end

[ss, cc] = meshgrid(sensory_infos, category_infos);

% Preallocate return values.
corrects = nan(size(ss));
slopes = nan(size(ss));
slopeErrors = nan(size(ss));

for i=1:numel(ss)
    params_copy = Model.setCategorySensoryInfo(params, cc(i), ss(i));
    
    % Construct UIDs for saving/loading with disk
    results_uid = Model.getModelStringID(params_copy);
    expfit_uid = ['PK-expfit-' optim_prefix '_' results_uid];

    % Set optimizable parameters to their optimal value
    if ~isempty(optim)
        for iVar=1:length(optim)
            params_copy.(optim{iVar}) = optim_params(i).(optim{iVar});
        end
    end
    
    % TODO - smarter setting of seed?
    params_copy.seed = randi(1000000000);
    
    % Run the model if needed
    [expFit, expErrors, results, correct_categories] = LoadOrRun(@runAndGetFit, {params_copy}, ...
        fullfile(memodir, expfit_uid));
    
    corrects(i) = mean(results.choices == correct_categories);
    slopes(i) = expFit(2);
    slopeErrors(i) = expErrors(2);
end

slopeErrors(isnan(slopeErrors)) = inf;
smoothSlopes = smoothn(slopes, 1./slopeErrors.^2);

colors = arrayfun(@(b) Model.betacolor(b, rb_range(1), rb_range(2)), linspace(rb_range(1), rb_range(2)), ...
    'UniformOutput', false);
cmap = vertcat(colors{:});

% gray_alphas = 1 - exp(-slopeErrors.^2/6);
gray_alphas = 1 - (1 + exp(-(corrects-.55)*100)).^-1;
gray_lum = 1;
gray_rgbs = gray_lum * ones([size(ss) 3]);

% Plot slope
fig = figure(); hold on;

% Slopes contourf
% contourf(smoothn(slopes), linspace(rb_range(1), rb_range(2), 11));
% axis image;
% colormap(cmap);

% Slopes image
imagesc('XData', sensory_infos, 'YData', category_infos, 'CData', smoothSlopes, rb_range);
axis image;
set(gca, 'YDir', 'normal');
colormap(cmap);
colorbar;

if exist('contour_pc', 'var')
    hold on;
    [~, threshold_pts] = Model.getThresholdPoints(category_infos, params, contour_pc, 50);
    plot(threshold_pts(:,1), threshold_pts(:,2), '-w', 'LineWidth', 2);

    if exist('contour_through_sc', 'var') && ~isempty(contour_through_sc)
        slopeVals = arrayfun(@(i) interp2(ss, cc, smoothSlopes, contour_through_sc(i,1), contour_through_sc(i,2)), ...
            1:size(contour_through_sc, 1));
        contour(ss, cc, smoothSlopes, slopeVals, 'k', 'LineWidth', 2);
    end
end

% Gray overlay
% h = image(gray_rgbs);
% set(h, 'AlphaData', gray_alphas);
% Labels n such

% Axis labels etc
set(gca, 'YTick', category_infos(round(linspace(1, length(category_infos), 5))));
set(gca, 'XTick', sensory_infos(round(linspace(1, length(sensory_infos), 5))));
xlabel('Sensory Info');
ylabel('Category Info');
title('C-S Space: PK Slope (\beta)');
figname = ['CSSlopes_' Model.getModelStringID(params, true) '.fig'];
saveas(fig, fullfile(savedir, figname));
end

function [expfit, experrors, runResults, correct_categories] = runAndGetFit(params)
uid = Model.getModelStringID(params);
runResults = LoadOrRun(@Model.runVectorized, {params}, fullfile(params.save_dir, uid));
[data, correct_categories] = Model.genDataWithParams(runResults.params);
% Fit PK to simulation (redo simulation on error)
err = true;
while err
    try
        [expfit, ~, experrors] = CustomRegression.ExponentialPK(data, runResults.choices == +1, 1);
        err = false;
    catch
        disp('trying again...');
        runResults = LoadOrRun(@Model.runVectorized, {params}, fullfile(params.save_dir, uid), ...
            '-recompute');
        [data, correct_categories] = Model.genDataWithParams(runResults.params);
    end
end
end