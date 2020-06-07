function visSubjectFitWrapper(subjectId, phaseStr, models, datadir, memodir)
if nargin < 3, models = {'is', 'itb', 'itb-gamma', 'ideal'}; end
if nargin < 4, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 5, memodir = fullfile(datadir, '..', 'Precomputed'); end
if ~iscell(models), models = {models}; end

[~, sigs, choices, ~, SubjectData] = GetSubjectDataForFitting(subjectId, 0.16);
switch lower(phaseStr)
    case 'hslc'
        phz = 1;
    case 'lshc'
        phz = 2;
    case 'both'
        phz = [1 2];
end
sigs = sigs(phz);
choices = choices(phz);
SubjectData = SubjectData(phz);

%% Load fit
for iMod=length(models):-1:1
    fits{iMod} = GetFit(subjectId, phaseStr, models{iMod}, true, datadir, memodir);
    complete(iMod) = ~isempty(fits{iMod});
    if complete(iMod), fields{iMod} = fits{iMod}(1).fit_fields; end
end

models = models(complete);
fits = fits(complete);
fields = fields(complete);

%% Simulate bestfit model forward and plot psychometric curve
figure;
for iMod=1:length(models)
    for iphz=1:length(phz)
        subplot(length(phz), length(models), (iphz-1)*length(models)+iMod); hold on;
        
        thisData = SubjectData{iphz};
        
        % Bin psychometric data for model and subject
        results = Model.runVectorized(fits{iMod}, sigs{iphz}/fits{iMod}.signal_scale);
        [xData, ~, bins] = unique(thisData.(thisData.pm_var));
        for iBin=length(xData):-1:1
            [pData(iBin,1), eData(iBin,:)] = binofit(sum(choices{iphz}(bins==iBin)==1), sum(bins==iBin));
            [modelPData(iBin,1), modelEData(iBin,:)] = binofit(sum(results.choices(bins==iBin)==1), sum(bins==iBin));
        end
        
        % Plot subject PM curve
        errorbar(xData, pData, pData-eData(:,1), eData(:,2)-pData, '-k');
        % Plot model PM curve
        errorbar(xData, modelPData, modelPData-modelEData(:,1), modelEData(:,2)-modelPData, '-r');
        
        if iphz*iMod==1, legend({'data', 'fit'}, 'location', 'best'); end
        title([models{iMod} ' :: [' num2str(phz(iphz)) ']']);
    end
end
sgtitle(subjectId);

%% Plot matrix of argmaxes
for iMod=1:length(models)
    figure;
    maxes = cellfun(@(fit) Fitting.getParamsFields(fit, fields{iMod}), fits{iMod}, 'uniformoutput', false);
    maxes = vertcat(maxes{:});
    ll = cellfun(@(fit) fit(1).ll, fits{iMod});
    ll_var = cellfun(@(fit) fit(1).ll_var, fits{iMod});
    [~, ibest] = max(ll);
    nF = length(fields{iMod});
    for iF=1:nF
        for jF=1:nF
            subplot(nF, nF, (iF-1)*nF+jF); hold on;
            if iF==jF
                histogram(maxes(:,iF), 50);
                plot(maxes(ibest,iF)*[1 1], ylim, '-k');
            else
                scatter(maxes(:,jF), maxes(:,iF), 15, ll, 'filled');
                plot(maxes(ibest,jF)*[1 1], ylim, '-k');
                plot(xlim, maxes(ibest,iF)*[1 1], '-k');
            end
            if iF==nF, xlabel(fields{iMod}{jF}); end
            if jF==1, ylabel(fields{iMod}{iF}); end
        end
    end
    sgtitle({models{iMod}, ...
        sprintf('Max LL = %.1f+/-%.1f', ll(ibest), sqrt(ll_var(ibest))), ...
        sprintf('Spread LL = %.1f+/-%.1f', mean(ll), std(ll))});
end

%% Error bar series of LL at each maximum
figure; hold on;
for iMod=1:length(models)
    ll = cellfun(@(fit) fit(1).ll, fits{iMod});
    ll_var = cellfun(@(fit) fit(1).ll_var, fits{iMod});
    [~,isrt] = sort(ll);
    errorbar(ll(isrt), sqrt(ll_var(isrt)));
    title(models{iMod});
    xlabel('run');
    ylabel('LL');
end
grid on;
legend(models, 'location', 'best');
end