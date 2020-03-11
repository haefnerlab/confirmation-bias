function [cs_fig, pk_fig] = plotCSPK(category_infos, sensory_infos, params, pk_hprs, pk_colormap, beta_range, clickpts)

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

% Keep in sync with plotCSSpace
figname = ['CSSpace_' Model.getModelStringID(params, true) '.fig'];

if ~exist(fullfile(savedir, figname), 'file')
    [~, cs_fig] = Model.plotCategorySensorySpace(category_infos, sensory_infos, params);
else
    cs_fig = openfig(fullfile(savedir, figname));
end

img_ax = findobj(cs_fig, 'Type', 'axes');
hold(img_ax, 'on');

if ~exist('clickpts', 'var')
    [x, y] = getpts(img_ax);
    npts = length(x);
    
    sens_pts = round(100*interp1(1:length(sensory_infos), sensory_infos, x, 'linear', 'extrap'))/100;
    cat_pts = round(100*interp1(1:length(category_infos), category_infos, y, 'linear', 'extrap'))/100;
    
    % 'snap' selected points to the given 'category_info' and 'sensory_info' grids.
    sens_pts = arrayfun(@(l) closest(sensory_infos, l), sens_pts);
    cat_pts = arrayfun(@(l) closest(category_infos, l), cat_pts);
else
    sens_pts = clickpts(:, 1);
    cat_pts = clickpts(:, 2);
    npts = length(sens_pts);
end

pk_fig = figure;
pk_ax = axes(pk_fig);

if nargin < 5 || isempty(pk_colormap)
    % create pk_colormap that is dark blue -> dark red
    fade = linspace(0, 1, npts)';
    colors = [fade*170/255, zeros(size(fade)), (1-fade)*170/255];
elseif isequal(pk_colormap, 'beta')
    colors = zeros(npts, 3);
else
    colors = pk_colormap(npts);
end

for i=1:length(sens_pts)
    s = sens_pts(i);
    c = cat_pts(i);
    
    params.category_info = c;
    params.sensory_info = s;
    params.p_match = c;
    params.var_s = Model.getEvidenceVariance(s);
    
    if isfield(params, 'gamma_min') && ~isempty(params.gamma_min)
        params.gamma = params.gamma_min + (params.gamma_max - params.gamma_min) * (1-c) * 2;
    end

    [weights, errors, para] = Model.plotPK(params, pk_hprs, {}, [], false);
    
    if isequal(pk_colormap, 'beta')
        % Compute color based on weights
        if ischar(pk_hprs) && startsWith(pk_hprs, 'exp')
            % para is abb, we want beta
            slope_param = para(2);
        elseif ischar(pk_hprs) && startsWith(pk_hprs, 'lin')
            % para is sob, we want slope
            slope_param = para(1);
        else
            % Auto-coloring enabled, but all we have to go off of is actual weights. Do post-hoc exp fit.
            ab = CustomRegression.expFit(weights(1:end-1), errors(1:end-1));
            slope_param = ab(2);
        end
        colors(i, :) = Model.betacolor(slope_param, beta_range(1), beta_range(2));
    end
    
    errors = errors / mean(weights(1:end-1));
    weights = weights / mean(weights(1:end-1));
    
    idx_s = length(sensory_infos) * (s-min(sensory_infos)) / (max(sensory_infos) - min(sensory_infos));
    idx_c = length(category_infos) * (c-min(category_infos)) / (max(category_infos) - min(category_infos));
    scatter(img_ax, idx_s, idx_c, 30, colors(i,:), 'filled');
    hold(pk_ax, 'on');
    errorbar(pk_ax, 1:params.frames, weights(1:end-1), errors(1:end-1), 'Color', colors(i,:), 'LineWidth', 2);
end

yl = ylim;
if yl(1) > 0
    ylim(pk_ax, [0 inf]);
end

%% Plot formatting

xlabel(img_ax, []);
ylabel(img_ax, []);
title(img_ax, []);
set(cs_fig, 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4]);
saveas(cs_fig, 'tmp_cspk_pts.fig');
saveas(cs_fig, 'tmp_cspk_pts.pdf');

set(pk_ax, 'XTick', [], 'YTick', [0 1], 'XAxisLocation', 'origin');
set(pk_fig, 'PaperSize', [6 4], 'PaperPosition', [0 0 6 4]);
saveas(pk_fig, 'tmp_cspk_lines.fig');
saveas(pk_fig, 'tmp_cspk_lines.pdf');

end