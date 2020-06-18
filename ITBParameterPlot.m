function [box_fig, hist_fig] = ITBParameterPlot(subjectIds, phases, plot_fields, log_flag, datadir, memodir)
if nargin < 3, plot_fields = {'prior_C', 'lapse', 'temperature', 'gamma', 'bound', 'noise'}; end
if nargin < 4, log_flag = false(size(plot_fields)); end
if nargin < 5, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 6, memodir = fullfile(datadir, '..', 'Precomputed'); end

% NOTE: doesn't work across phases w/ different fields!


allsamples = zeros(0, length(plot_fields));
allgroups = [];
for iSub=1:length(subjectIds)
    for iPhz=1:length(phases)
        [chains, fields{iPhz}, loglikes{iSub,iPhz}, ~, ~, ~, ~] = LoadOrRun(@GetITBPosteriorSamples, ...
            {subjectIds{iSub}, phases{iPhz}, 0, true, 1:12, datadir, memodir}, ...
            fullfile('tmp-mh', ['ITB-cache-' subjectIds{iSub} '-' phases{iPhz} '.mat']));

        samples = vertcat(chains{:});
        values{iSub,iPhz} = zeros(size(samples,1), length(plot_fields));
        for iPara=length(plot_fields):-1:1
            iF = cellfun(@(f) contains(f, plot_fields{iPara}, 'ignorecase', true), fields{iPhz});
            if log_flag(iPara)
                values{iSub,iPhz}(:,iPara) = log10(samples(:,iF));
            else
                values{iSub,iPhz}(:,iPara) = samples(:,iF);
            end
        end
    end
end

%% Get chain-reweighted moments and quantiles
range_lo = min(vertcat(values{:}));
range_hi = max(vertcat(values{:}));
for iPara=length(plot_fields):-1:1
    edges(iPara, :) = linspace(range_lo(iPara), range_hi(iPara), 101);
end

for iSub=length(subjectIds):-1:1
    for iPhz=length(phases):-1:1
        [~, quantiles{iSub,iPhz}, pmfs{iSub,iPhz}] = ...
            sampleQuantilesReweightChains(values{iSub,iPhz}, loglikes{iSub,iPhz}, [.025 .25 .5 .75 .975], edges);
    end
end

%% Scatters
scatter_fig = figure;
for iPara=1:length(plot_fields)
    subplot(1, length(plot_fields), iPara); hold on;
    for iSub=1:length(subjectIds)
        [m1, l1, u1] = deal(quantiles{iSub,1}(iPara,3), quantiles{iSub,1}(iPara,1), quantiles{iSub,1}(iPara,5));
        [m2, l2, u2] = deal(quantiles{iSub,2}(iPara,3), quantiles{iSub,2}(iPara,1), quantiles{iSub,2}(iPara,5));
        errorbar(m1, m2, m2-l2, u2-m2, m1-l1, u1-m1, '.b');
    end
    axis equal; axis square; grid on;
    title(plot_fields{iPara});
    xlabel(phases{1}); ylabel(phases{2});
    uistack(plot(xlim, xlim, '-k', 'HandleVisibility', 'off'), 'bottom');
end

%% Marginal histograms figure
hist_fig = figure;
for iSub=1:length(subjectIds)
    for iPara=1:length(plot_fields)
        % One plot per subject per parameter. Stacked vertically across subjects so each parameter
        % shares an x-axis. Overlay histograms per phase
        subplot(length(subjectIds), length(plot_fields), (iSub-1)*length(plot_fields)+iPara);
        hold on;
        
        ctrs = (edges(iPara,1:end-1) + edges(iPara,2:end))/2;
        bar(ctrs, pmfs{iSub, 1}(iPara, :), 1, 'FaceAlpha', .5, 'EdgeAlpha', 0);
        bar(ctrs, pmfs{iSub, 2}(iPara, :), 1, 'FaceAlpha', .5, 'EdgeAlpha', 0);
        
        set(gca, 'YTick', 0);
        if iSub < length(subjectIds)
            set(gca, 'XTickLabel', []);
        else
            xlabel(plot_fields{iPara});
        end
        if iPara == 1
            ylabel(shortname(subjectIds{iSub}));
        end
    end
end

%% Box plot figure
box_fig = figure;
for iPara=1:length(plot_fields)
    subplot(1, length(plot_fields), iPara); hold on;

    % 95% CI as a thin line
    for iSub=1:length(subjectIds)
        plot(quantiles{iSub,1}(iPara,[1 5]), iSub*[1 1], '-', 'Color', 'r', 'LineWidth', 1);
        plot(quantiles{iSub,2}(iPara,[1 5]), iSub*[1 1]+.25, '-', 'Color', 'b', 'LineWidth', 1);
    end

    % 50% CI as a thick line
    for iSub=1:length(subjectIds)
        plot(quantiles{iSub,1}(iPara,[2 4]), iSub*[1 1], '-', 'Color', 'r', 'LineWidth', 2);
        plot(quantiles{iSub,2}(iPara,[2 4]), iSub*[1 1]+.25, '-', 'Color', 'b', 'LineWidth', 2);
    end

    % Median as a black dot
    for iSub=1:length(subjectIds)
        plot(quantiles{iSub,1}(iPara,3), iSub, '.', 'Color', 'k');
        plot(quantiles{iSub,2}(iPara,3), iSub+.25, '.', 'Color', 'k');
    end

    title(plot_fields{iPara});
    set(gca, 'YTick', 1:length(subjectIds), 'YTickLabel', cellfun(@shortname, subjectIds, 'uniformoutput', false));
    grid on;
end

end