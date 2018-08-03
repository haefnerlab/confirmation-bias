function plotCSBias(category_infos, sensory_infos, params, clickpts)
savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

% Keep in sync with plotCSSpace
figname = ['CSSpace_' Model.getModelStringID(params, true) '.fig'];

if ~exist(fullfile(savedir, figname), 'file')
    Model.plotCategorySensorySpace(category_infos, sensory_infos, params);
    close all;
end

cs_fig = openfig(fullfile(savedir, figname));
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

% create pk_colormap that is dark blue -> dark red
fade = linspace(0, 1, npts)';
colors = [fade*170/255, zeros(size(fade)), (1-fade)*170/255];

bias_fig = figure;

for ipt=1:npts
    s = sens_pts(ipt);
    c = cat_pts(ipt);
    [~, il] = min(abs(sensory_infos - s));
    [~, ip] = min(abs(category_infos - c));
    
    params.category_info = c;
    params.sensory_info = s;
    params.p_match = c;
    params.var_s = Model.getEvidenceVariance(s);
    
    figure(bias_fig);
    subplotsquare(npts, ipt);
    Model.plotTrueBiasCorrection(params, {'Color', colors(ipt, :)});
    
    scatter(img_ax, il, ip, 50, colors(ipt,:), 'filled');
end

end