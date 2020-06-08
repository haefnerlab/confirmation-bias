function ITBParameterPlot(subjectIds, phases, params, datadir, memodir)
if nargin < 3, params = {'prior_C', 'temperature', 'gamma', 'bound', 'noise'}; end
if nargin < 4, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 5, memodir = fullfile(datadir, '..', 'Precomputed'); end

% NOTE: doesn't work across phases w/ different fields!

allsamples = [];
allgroups = [];
for iSub=1:length(subjectIds)
    for iPhz=1:length(phases)
        [samples, fields] = GetITBPosteriorSamples(subjectIds{iSub}, phases{iPhz}, 0, true, datadir, memodir);
        
        groupId = (iSub-1)*length(phases)+iPhz;
        allsamples = vertcat(allsamples, samples);
        allgroups = vertcat(allgroups, groupId*ones(size(samples,1),1));
        groupname{groupId} = sprintf('%s [%s]', subjectIds{iSub}, phases{iPhz});
    end
end
        
for iPara=1:length(params)
    subplot(1,length(params), iPara); hold on;
    iF = cellfun(@(f) contains(f, params{iPara}, 'ignorecase', true), fields);
    if any(iF)
        boxplot(allsamples(:,iF), allgroups, 'plotstyle', 'compact', 'symbol', '', 'colorgroup', mod(allgroups, length(phases)));
    end
    title(params{iPara});
    set(gca, 'XTick', 1:length(groupname), 'XTickLabel', groupname);
    xtickangle(90);
    grid on;
end

end