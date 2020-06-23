function [scatter_fig, hist_fig, box_fig] = ITBParameterPlot(subjectIds, phases, plot_fields, log_flag, colors, markers, datadir, memodir)
if nargin < 3, plot_fields = {'prior_C', 'lapse', 'temperature', 'gamma', 'bound', 'noise'}; end
if nargin < 4, log_flag = false(size(plot_fields)); end
if nargin < 5, colors = lines(length(plot_fields)); end
if nargin < 6, markers = repmat('o', 1, length(plot_fields)); end
if nargin < 7, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 8, memodir = fullfile(datadir, '..', 'Precomputed'); end

for iSub=1:length(subjectIds)
    for iPhz=1:length(phases)
        [chains, fields{iPhz}, ~, logpost{iSub,iPhz}, ~, ~, ~, ~] = LoadOrRun(@GetITBPosteriorSamples, ...
            {subjectIds{iSub}, phases{iPhz}, 0, true, 1:12, datadir, memodir}, ...
            fullfile('tmp-mh', ['ITB-cache-' subjectIds{iSub} '-' phases{iPhz} '.mat']));
        
        % Do thinning down to 500 samples per chain
        samples = vertcat(chains{:});
        istart = 1;
        thin_idx = cell(size(chains));
        for iChain=1:length(chains)
            iend = istart + length(logpost{iSub,iPhz}{iChain}) - 1;
            assert(length(logpost{iSub,iPhz}{iChain}) >= 500);
            thin_idx{iChain} = round(linspace(istart, iend, 500));
            istart = iend+1;
        end
        thin_idx = horzcat(thin_idx{:});
        thin_samples = samples(thin_idx, :);

        % Convert from [S x #params] all-samples to [S x length(plot_fields)] values we care about
        values{iSub,iPhz} = zeros(size(thin_samples,1), length(plot_fields));
        for iPara=length(plot_fields):-1:1
            iF = cellfun(@(f) contains(f, plot_fields{iPara}, 'ignorecase', true), fields{iPhz});
            if any(iF)
                values{iSub,iPhz}(:,iPara) = thin_samples(:,iF);
            else
                values{iSub,iPhz}(:,iPara) = nan(size(thin_samples,1),1);
            end
        end
    end
end

dists = Fitting.defaultDistributions(plot_fields);
lb = cellfun(@(f) dists.(f).plb, plot_fields);
ub = cellfun(@(f) dists.(f).pub, plot_fields);

%% Get chain-reweighted moments and quantiles
range_lo = min(vertcat(values{:}));
range_hi = max(vertcat(values{:}));
for iPara=length(plot_fields):-1:1
    if log_flag(iPara)
        edges(iPara, :) = logspace(log10(range_lo(iPara)), log10(range_hi(iPara)), 101);
    else
        edges(iPara, :) = linspace(range_lo(iPara), range_hi(iPara), 101);
    end
end

% [1 7] are 95% interval, [2 6] are 68% interval. [3 5] are 50% interval, and [4] is median
quants = [.025 .1587 .25 .5 .75 .8414 .975];
for iSub=length(subjectIds):-1:1
    for iPhz=length(phases):-1:1
        % EXPERIMENTAL: reweight chains by average density. Call helper to get adjusted quantiles and PMFs
        % [~, quantiles{iSub,iPhz}, pmfs{iSub,iPhz}] = ...
        %     sampleQuantilesReweightChains(values{iSub,iPhz}, logpost{iSub,iPhz}, quants, edges);
        
        % STANDARD: get quantiles and PMFs the old fashioned way
        for iQ=length(quants):-1:1
            quantiles{iSub,iPhz}(:,iQ) = quantile(values{iSub,iPhz}, quants(iQ));
        end

        for iPara=length(plot_fields):-1:1
            pmfs{iSub,iPhz}(iPara,:) = histcounts(values{iSub,iPhz}(:,iPara), edges(iPara,:));
        end
    end
end

%% Scatters @ 68% confidence intervals
scatter_fig = figure;
for iPara=1:length(plot_fields)
    subplot(1, length(plot_fields), iPara); hold on;
    for iSub=1:length(subjectIds)
        [m1, l1, u1] = deal(quantiles{iSub,1}(iPara,4), quantiles{iSub,1}(iPara,2), quantiles{iSub,1}(iPara,6));
        [m2, l2, u2] = deal(quantiles{iSub,2}(iPara,5), quantiles{iSub,2}(iPara,2), quantiles{iSub,2}(iPara,6));
        errorbar(m1, m2, m2-l2, u2-m2, m1-l1, u1-m1, markers(iPara), 'Color', colors(iPara, :), 'CapSize', 0, 'MarkerFaceColor', colors(iPara,:), 'MarkerEdgeColor', [1 1 1]);
        if log_flag(iPara)
            set(gca, 'XScale', 'log', 'YScale', 'log');
        end
    end
    axis equal; axis square; %grid on;
    title(plot_fields{iPara});
    xlabel(phases{1}); ylabel(phases{2});
    uistack(plot([lb(iPara) ub(iPara)], [lb(iPara) ub(iPara)], '-k', 'HandleVisibility', 'off'), 'bottom');
    xlim([lb(iPara) ub(iPara)]); ylim([lb(iPara) ub(iPara)]);
end

%% Marginal histograms figure
hist_fig = figure;
for iSub=1:length(subjectIds)
    for iPara=1:length(plot_fields)
        % One plot per subject per parameter. Stacked vertically across subjects so each parameter
        % shares an x-axis. Overlay histograms per phase
        subplot(length(subjectIds), length(plot_fields), (iSub-1)*length(plot_fields)+iPara);
        hold on;
        
        if log_flag(iPara)
            ctrs = (log10(edges(iPara,1:end-1)) + log10(edges(iPara,2:end)))/2;
            xticks = log10(linspace(10.^ctrs(1), 10.^ctrs(end), 20));
            xlab = 10.^xticks;
        else
            ctrs = (edges(iPara,1:end-1) + edges(iPara,2:end))/2;
            xticks = linspace(ctrs(1), ctrs(end), 20);
            xlab = xticks;
        end
        bar(ctrs, pmfs{iSub, 1}(iPara, :), 1, 'FaceAlpha', .5, 'EdgeAlpha', 0);
        bar(ctrs, pmfs{iSub, 2}(iPara, :), 1, 'FaceAlpha', .5, 'EdgeAlpha', 0);
        
        set(gca, 'YTick', 0, 'XTick', xticks);
        if iSub < length(subjectIds)
            set(gca, 'XTickLabel', []);
        else
            set(gca, 'XTickLabel', arrayfun(@(n) num2str(n,2), xlab, 'uniformoutput', false));
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
        plot(quantiles{iSub,1}(iPara,[1 7]), iSub*[1 1], '-', 'Color', 'r', 'LineWidth', 1);
        plot(quantiles{iSub,2}(iPara,[1 7]), iSub*[1 1]+.25, '-', 'Color', 'b', 'LineWidth', 1);
    end
    
    % 50% CI as a thick line
    for iSub=1:length(subjectIds)
        plot(quantiles{iSub,1}(iPara,[3 5]), iSub*[1 1], '-', 'Color', 'r', 'LineWidth', 2);
        plot(quantiles{iSub,2}(iPara,[3 5]), iSub*[1 1]+.25, '-', 'Color', 'b', 'LineWidth', 2);
    end
    
    % Median as a black dot
    for iSub=1:length(subjectIds)
        plot(quantiles{iSub,1}(iPara,4), iSub, '.', 'Color', 'k');
        plot(quantiles{iSub,2}(iPara,4), iSub+.25, '.', 'Color', 'k');
    end
    
    if log_flag(iPara), set(gca, 'XScale', 'log'); end
    
    title(plot_fields{iPara});
    set(gca, 'YTick', 1:length(subjectIds), 'YTickLabel', cellfun(@shortname, subjectIds, 'uniformoutput', false));
    if log_flag(iPara)
        xl = xlim;
        set(gca, 'XMinorGrid', 'off', 'xtick', 10.^(ceil(log10(xl(1))):floor(log10(xl(2)))))
    end
    ylim([0.75 length(subjectIds)+0.5]);
    grid on;
end

end