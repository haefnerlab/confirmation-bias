function plotCSPK(category_infos, sensory_infos, params, pk_hprs)

savedir = fullfile('+VariationalModel', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 4, pk_hprs = [1 0 10]; end

% Keep in sync with plotCSSpace

if params.noise > 0
    noise_str = ['_' num2str(params.noise, 2)];
else
    noise_str = '';
end

figname = sprintf('CSSpace_%dx%d_vx%.2f_pC%.2f_u%d_gam%.2f%s.fig', params.trials, params.frames, ...
    params.var_x, params.prior_C, params.updates, params.gamma, noise_str);

if ~exist(fullfile(savedir, figname), 'file')
    VariationalModel.plotCategorySensorySpace(category_infos, sensory_infos, params);
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
    params.var_s = SamplingModel.getEvidenceVariance(s);
    [weights, errors, tmp_fig] = VariationalModel.plotPK(params, pk_hprs);
    close(tmp_fig);
    hold(pk_ax, 'on');
    errorbar(1:params.frames, weights(1:end-1), errors(1:end-1), 'Color', colors(i,:), 'LineWidth', 2);
end

yl = ylim;
if yl(1) > 0
    ylim(pk_ax, [0 inf]);
end

end