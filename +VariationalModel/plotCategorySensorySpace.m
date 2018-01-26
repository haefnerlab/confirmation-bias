function correct = plotCategorySensorySpace(category_infos, sensory_infos, params)
%PLOTCATEGORYSENSORYSPACE make category_info vs sensory_info plots for the given params.

savedir = fullfile('+VariationalModel', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

[ss, cc] = meshgrid(sensory_infos, category_infos);

% Preallocate return values.
correct = nan(size(ss));

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
    [~, ~, choices] = LoadOrRun(@VariationalModel.runLatentZNormalX, {params_copy}, ...
        fullfile(params.save_dir, results_uid));
    correct(i) = mean(choices == +1);
end

% Plot percent correct
figure();
imagesc(correct, [0.5 1.0]); axis image; colorbar;
% Add contour line at threshold
hold on; contour(imfilter(correct, smoothkernel(9), 'replicate'), [0.7 0.7], '-w', 'LineWidth', 2);
% Add labels and axis ticks
category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
set(gca, 'YTick', category_tick_indices);
set(gca, 'XTick', sensory_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('Sensory Info');
ylabel('Category Info');
title('C-S Space: Percent Correct');

figname = sprintf('CSSpace_%dx%d_vx%.2f_pC%.2f_u%d_gam%.2f.fig', params.trials, params.frames, ...
    params.var_x, params.prior_C, params.updates, params.gamma);
saveas(gcf, fullfile(savedir, figname));
end

function kernel = smoothkernel(n)
x = linspace(-2,2,n);
[xx,yy] = meshgrid(x,x);
kernel = exp(-(xx.^2 + yy.^2));
kernel = kernel / sum(kernel(:));
end