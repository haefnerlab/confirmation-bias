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
        [chains, fields, ~, ~, ~, ~, ~] = LoadOrRun(@GetITBPosteriorSamples, ...
            {subjectIds{iSub}, phases{iPhz}, 0, true, 1:12, datadir, memodir}, ...
            fullfile('tmp-mh', ['ITB-cache-' subjectIds{iSub} '-' phases{iPhz} '.mat']));
        samples = vertcat(chains{:});
        
        % Format expected by 'boxplot' below is a massive concatenation of all samples, but with a
        % 'group' id per row indicating subject x condition.
        groupId = (iSub-1)*length(phases)+iPhz;
        allsamples(end+1:end+size(samples,1), 1:size(samples,2)) = samples;
        allgroups = vertcat(allgroups, groupId*ones(size(samples,1),1));
        groupname{groupId} = sprintf('%s [%s]', shortname(subjectIds{iSub}), phases{iPhz});
    end
end

%% Pull out just 'plot_fields' columns (and maybe log transform)
for iPara=length(plot_fields):-1:1
    iF = cellfun(@(f) contains(f, plot_fields{iPara}, 'ignorecase', true), fields);
    if log_flag(iPara)
        allvalues(:,iPara) = log10(allsamples(:,iF));
        plot_fields{iPara} = ['log_{10}(' plot_fields{iPara} ')'];
    else
        allvalues(:,iPara) = allsamples(:,iF);
    end
end

range_lo = quantile(allvalues, .01);
range_hi = quantile(allvalues, .99);

%% Scatters
scatter_fig = figure;
for iPara=1:length(plot_fields)
end

%% Marginal histograms figure
hist_fig = figure;
for iSub=1:length(subjectIds)
    for iPara=1:length(plot_fields)
        % One plot per subject per parameter. Stacked vertically across subjects so each parameter
        % shares an x-axis. Overlay histograms per phase
        subplot(length(subjectIds), length(plot_fields), (iSub-1)*length(plot_fields)+iPara);
        hold on;
        
        bin_edges = linspace(range_lo(iPara), range_hi(iPara), 101);
        for iPhz=1:length(phases)
            % Format expected by 'boxplot' below is a massive concatenation of all samples, but with a
            % 'group' id per row indicating subject x condition.
            groupId = (iSub-1)*length(phases)+iPhz;
            
            histogram(allvalues(allgroups == groupId, iPara), bin_edges, 'EdgeAlpha', 0.1);
        end
        
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
    boxplot(allvalues(:,iPara), allgroups, 'plotstyle', 'compact', 'symbol', '', 'colorgroup', mod(allgroups, length(phases)), 'orientation', 'horizontal');
    title(plot_fields{iPara});
    set(gca, 'YTick', 1:length(groupname), 'YTickLabel', groupname);
    grid on;
end

end