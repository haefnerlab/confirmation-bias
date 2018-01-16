function plotCSPK(category_infos, sensory_infos, params, ideal_observer, pk_hprs, optimize, optim_grid_size, pk_colormap)

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end
if nargin < 6, optimize = {}; end
if nargin < 7, optim_grid_size = 11; end

optim_prefix = Model.getOptimPrefix(optimize, optim_grid_size);

% Keep in sync with plotCSSpace
if ~ideal_observer
    figname = sprintf('CSSpace_%dx%d_%s_vx%.2f_pC%.2f_gam%.2f_ns%d_b%d_%d.fig', ...
        params.trials, params.frames, optim_prefix, params.var_x, params.prior_C, ...
        params.gamma, params.samples, params.batch, params.importance_norm);
else
    figname = sprintf('CSSpace_%dx%d_vx%.2f_ideal.fig', params.trials, params.frames, params.var_x);
end

if ~exist(fullfile(savedir, figname), 'file')
    Model.plotCategorySensorySpace(category_infos, sensory_infos, params, ideal_observer, optimize, optim_grid_size);
    close all;
end

cs_fig = openfig(fullfile(savedir, figname));
img_ax = findobj(cs_fig, 'Type', 'axes');
hold(img_ax, 'on');

[x, y] = getpts(img_ax);
npts = length(x);

sens_pts = round(100*interp1(1:length(sensory_infos), sensory_infos, x, 'linear', 'extrap'))/100;
cat_pts = round(100*interp1(1:length(category_infos), category_infos, y, 'linear', 'extrap'))/100;

% 'snap' selected points to the given 'category_info' and 'sensory_info' grids.
sens_pts = arrayfun(@(l) closest(sensory_infos, l), sens_pts);
cat_pts = arrayfun(@(l) closest(category_infos, l), cat_pts);

pk_fig = figure;
pk_ax = axes(pk_fig);

if nargin < 8
    % create pk_colormap that is dark blue -> dark red
    fade = linspace(0, 1, npts)';
    colors = [fade*170/255, zeros(size(fade)), (1-fade)*170/255];
else
    colors = pk_colormap(npts);
end

for i=1:length(sens_pts)
    s = sens_pts(i);
    c = cat_pts(i);
    [~, il] = min(abs(sensory_infos - s));
    [~, ip] = min(abs(category_infos - c));
    scatter(img_ax, il, ip, 50, colors(i,:), 'filled');
    
    params.category_info = c;
    params.sensory_info = s;
    params.p_match = c;
    params.var_e = Model.getEvidenceVariance(s);
    [weights, errors, tmp_fig] = Model.plotSamplingPK(params, pk_hprs, ideal_observer, optimize, optim_grid_size);
    close(tmp_fig);
    hold(pk_ax, 'on');
    errorbar(1:params.frames, weights(1:end-1), errors(1:end-1), 'Color', colors(i,:), 'LineWidth', 2);
end

yl = ylim;
if yl(1) > 0
    ylim(pk_ax, [0 inf]);
end

end