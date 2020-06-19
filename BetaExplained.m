function [fig_bar, fig_scatter] = BetaExplained(subjectIds, plotstyle, datadir, memodir)
if nargin < 2, plotstyle = 'stacked'; end
if nargin < 3, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 4, memodir = fullfile(datadir, '..', 'Precomputed'); end

% GBN = gamma/bound/noise. inv = inverted.
assert(ismember(lower(plotstyle), {'gbn', 'gbn-inv', 'separate'}));

% phase_names = {'hslc', 'lshc', 'both', 'both'};
% plot_titles = {'HSLC', 'LSHC', 'Both [HSLC]', 'Both [LSHC]'};
% data_idxs = [1 1 1 2];
phase_names = {'hslc', 'lshc'};
plot_titles = {'HSLC', 'LSHC'};
data_idxs = [1 1];
test_fields = {'gamma', 'bound', 'noise'};

%% Load posterior samples per subject per condition
for iSub=length(subjectIds):-1:1
    for iPhz=length(phase_names):-1:1
        % Load posterior samples over all parameters
        [chains, fields, loglike, ~, params, sigs, choices] = LoadOrRun(@GetITBPosteriorSamples, ...
            {subjectIds{iSub}, phase_names{iPhz}, 0, true, 1:12, datadir, memodir}, ...
            fullfile('tmp-mh', ['ITB-cache-' subjectIds{iSub} '-' phase_names{iPhz} '.mat']), '-verbose');
        samples = vertcat(chains{:});
        sigs = sigs{data_idxs(iPhz)};
        choices = choices{data_idxs(iPhz)};
        params = params(data_idxs(iPhz));
        if strcmp(phase_names{iPhz}, 'both')            
            % If 'both' condition, we have duplicate fields and params for HSLC and LSHC conditions.
            % Pare down fields and samples to just those relevant to the currently tested condition.
            this_condition = data_idxs(iPhz);
            other_condition = 3-this_condition; % 2 --> 1 or 1 --> 2
            fields_other_condition = cellfun(@(f) endsWith(f, sprintf('_%d', other_condition)), fields);
            
            samples = samples(:, ~fields_other_condition);
            fields = fields(~fields_other_condition);
            
            % Remove '_1' or '_2' suffix from fields for this condition 'as if' just fit to the
            % current condition.
            fields_this_condition = cellfun(@(f) endsWith(f, sprintf('_%d', this_condition)), fields);
            fields(fields_this_condition) = cellfun(@(f) f(1:end-2), fields(fields_this_condition), 'uniformoutput', false);
        end
        
        % Thin samples since they are generally autocorrelated and this saves on redundant
        % computation effort
        w = sampleQuantilesReweightChains(samples, loglike, []);
        thin_idx = round(linspace(1, size(samples,1), min(500, size(samples,1))));
        samples = samples(thin_idx, :);
        weights{iSub, iPhz} = w(thin_idx);
        
        % Do ablation test on 'test' parameters. Since this involves simulating but not
        % marginalizing over the model simulations, set model to 'itb' rather than 'itb-int' for
        % speed.
        assert(strcmp(params.model, 'itb-int'), 'sanity failed');
        params.model = 'itb';
        [betaExplained{iSub, iPhz}, ~, ablations{iPhz}] = LoadOrRun(@PKShapeExplained, ...
            {sigs, params, test_fields, fields, samples}, ...
            fullfile('tmp-mh', ['PK-beta-cache-' subjectIds{iSub} '-' phase_names{iPhz} '-' num2str(data_idxs(iPhz)) '-' strjoin(test_fields) '.mat']), '-verbose');
        
        % Get 'true' slope of exponential PK for this subject x phase
        memo_name = ['Boot-ExpPK-' subjectIds{iSub} '-' phase_names{iPhz} '-fitting.mat'];
        [~, ~, ~, ~, ~, abb{iSub,iPhz}] = LoadOrRun(@BootstrapExponentialWeightsGabor, {sigs, choices, 10000, true}, ...
            fullfile(memodir, memo_name));
    end
end

idx_full = cellfun(@isempty, ablations{iPhz});
idx_null = cellfun(@(abl) isequal(abl, test_fields), ablations{iPhz});
idx_g = cellfun(@(abl) isequal(abl, {'gamma'}), ablations{iPhz});
idx_bn = cellfun(@(abl) isequal(abl, {'bound', 'noise'}), ablations{iPhz});

%% Scatter plots

% Hack... this section also only works for specific hypothesis tests on g/b/n
assert(isequal(test_fields, {'gamma', 'bound', 'noise'}));

color0 = [0.1 0.1 0.1];
color1 = [0 .5 0];
color2 = [.4 0 .4];
                    
fig_scatter = figure;
for iPhz=1:length(phase_names)
    subplot(1,length(phase_names),iPhz);
    hold on;
    
    % populate true beta
    mtrue = zeros(size(subjectIds)); ltrue = zeros(size(subjectIds)); utrue = zeros(size(subjectIds));
    for iSub=1:length(subjectIds)
        [mtrue(iSub), ltrue(iSub), utrue(iSub)] = meanci(abb{iSub,iPhz}(:,2), 0.68);
    end
    % populate full model beta
    mfull = zeros(size(subjectIds)); lfull = zeros(size(subjectIds)); ufull = zeros(size(subjectIds));
    for iSub=1:length(subjectIds)
        [mfull(iSub), lfull(iSub), ufull(iSub)] = quantileWeight(betaExplained{iSub,iPhz}(:,idx_full), 0.68, weights{iSub,iPhz});
    end
    % populate gamma-ablated beta
    mgam = zeros(size(subjectIds)); lgam = zeros(size(subjectIds)); ugam = zeros(size(subjectIds));
    for iSub=1:length(subjectIds)
        [mgam(iSub), lgam(iSub), ugam(iSub)] = quantileWeight(betaExplained{iSub,iPhz}(:,idx_g), 0.68, weights{iSub,iPhz});
    end
    % populate bn-ablated beta
    mbn = zeros(size(subjectIds)); lbn = zeros(size(subjectIds)); ubn = zeros(size(subjectIds));
    for iSub=1:length(subjectIds)
        [mbn(iSub), lbn(iSub), ubn(iSub)] = quantileWeight(betaExplained{iSub,iPhz}(:,idx_bn), 0.68, weights{iSub,iPhz});
    end
    
    % Plot full model vs true beta 
    errorbar(mtrue, mfull, mfull-lfull, ufull-mfull, mtrue-ltrue, utrue-mtrue, 'o', 'Color', color0, 'MarkerFaceColor', color0, 'MarkerEdgeColor', [1 1 1]);
    % Plot gamma-ablated model vs true beta
    errorbar(mtrue, mgam, mgam-lgam, ugam-mgam, mtrue-ltrue, utrue-mtrue, '^', 'Color', color1, 'MarkerFaceColor', color1, 'MarkerEdgeColor', [1 1 1]);
    % Plot bn-ablated model vs true beta
    errorbar(mtrue, mbn, mbn-lbn, ubn-mbn, mtrue-ltrue, utrue-mtrue, 's', 'Color', color2, 'MarkerFaceColor', color2, 'MarkerEdgeColor', [1 1 1]);
    
    xlim([-.5 .5]); ylim([-.5 .5]);
    axis square;
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'XTick', -.5:.1:.5, 'YTick', -.5:.1:.5);
    grid on;
    uistack(plot([-.5 .5], [-.5 .5], '-k', 'HandleVisibility', 'off'), 'bottom');
    legend({'full', 'leak (\gamma) = 0', 'bound = \infty'}, 'location', 'best');
            
    title(plot_titles{iPhz});
end
drawnow;
pause(0.1);

% %% Decomposition patch plots
% 
% % Hack... this section also only works for specific hypothesis tests on g/b/n
% assert(isequal(test_fields, {'gamma', 'bound', 'noise'}));
% 
% color0 = [0.1 0.1 0.1];
% color1 = [0 .5 0];
% color2 = [.4 0 .4];
%                     
% fig_decompose = figure;
% for iPhz=1:length(phase_names)
%     subplot(1,length(phase_names),iPhz);
%     hold on;
%     
%     % get true, full, g, bn, and sum betas
%     mtrue = zeros(size(subjectIds)); ltrue = zeros(size(subjectIds)); utrue = zeros(size(subjectIds));
%     for iSub=1:length(subjectIds)
%         [mtrue(iSub), ltrue(iSub), utrue(iSub)] = meanci(abb{iSub,iPhz}(:,2), 0.68);
% 
%         gvals = betaExplained{iSub,iPhz}(:,idx_g);
%         bvals = betaExplained{iSub,iPhz}(:,idx_bn);
%         sumvals = gvals + bvals;
%         [sumvals, isrt] = sort(sumvals);
%         
% %         rescale = fvals ./ sumvals;
% %         gvals = gvals .* rescale;
% %         bvals = bvals .* rescale;
%         
%         gpatchx = [linspace(-.5, .5, length(gvals)), .5, -.5, -.5];
%         gpatchy = [gvals(isrt)' 0 0 gvals(1)];
%         bpatchx = [linspace(-.5, .5, length(gvals)) linspace(.5, -.5, length(bvals))];
%         bpatchy = [gvals(isrt)' fliplr(sumvals')];
%         patch(iSub+gpatchx,gpatchy,color2,'FaceAlpha',.5);
%         patch(iSub+bpatchx,bpatchy,color1,'FaceAlpha',.5);
%     end
%             
%     title(plot_titles{iPhz});
% end
% drawnow;
% pause(0.1);

%% Bar plots

fig_bar = figure;
for iSub=1:length(subjectIds)
    for iPhz=1:length(phase_names)
        if ~isempty(betaExplained{iSub, iPhz})
            subplot(length(subjectIds), length(phase_names), (iSub-1)*length(phase_names) + iPhz);
            hold on;
            switch lower(plotstyle)
                case 'separate'
                    [meanB, loB, hiB] = quantileWeight(betaExplained{iSub, iPhz}, .68, weights{iSub,iPhz});
                    bar(1:length(meanB), meanB, 'FaceColor', [.9 .9 .9]);
                    errorbar(1:length(meanB), meanB, meanB-loB, hiB-meanB, 'ok');
                    
                    [m,l,u] = meanci(abb{iSub,iPhz}(:,2), 0.68);
                    bar(length(meanB)+1, m, 'FaceColor', [1 1 1]);
                    errorbar(length(meanB)+1, m, m-l, u-m, 'ok');
                    
                    if inverted, warning('not implemented'); end
                    
                    abl = ablations{iPhz};
                    names = cellfun(@(a) strjoin(a, '+'), abl, 'uniformoutput', false);
                    names = [names {'true'}];
                    set(gca, 'XTick', 1:length(names), 'XTickLabel', names);
                    xtickangle(60);
                case 'gbn'
                    % In non-inverted 'normal' case, contribution is beta added relative to
                    % all-ablated

                    % This is a targeted analysis comparing gamma/bound/noise... use 'separate' flag
                    % for more generic parameter combinations
                    assert(isequal(test_fields, {'gamma', 'bound', 'noise'}));
                                        
                    % Get betas for 'full model', gamma-ablated, and bn-ablated
                    [barData(1), errData(1,1), errData(2,1)] = quantileWeight(betaExplained{iSub,iPhz}(:,idx_g), 0.68, weights{iSub,iPhz});
                    [barData(2), errData(1,2), errData(2,2)] = quantileWeight(betaExplained{iSub,iPhz}(:,idx_bn), 0.68, weights{iSub,iPhz});
                    [barData(3), errData(1,3), errData(2,3)] = quantileWeight(betaExplained{iSub,iPhz}(:,idx_full), 0.68, weights{iSub,iPhz});
                    [barData(4), errData(1,4), errData(2,4)] = meanci(abb{iSub,iPhz}(:,2), 0.68);
                    
                    h = bar([1 nan], [barData; nan(1,4)]);
                    drawnow; % populate h.xoffset
                    errorbar(h(1).XData(1) + [h.XOffset], barData, errData(1,:), errData(2,:), '.k');
                    h(1).FaceColor = color1;
                    h(2).FaceColor = color2;
                    h(3).FaceColor = color0;
                    h(4).FaceColor = [1 1 1];
                    
                    if iSub*iPhz == 1, legend({'leak (\gamma) = 0', 'bound = \infty', 'full', 'true'}); end
                case 'gbn-inv'
                    % In 'inverted' case, contribution of each parameter is the beta that is
                    % lost when that parameter is ablated

                    % This is a targeted analysis comparing gamma/bound/noise... use 'separate' flag
                    % for more generic parameter combinations
                    assert(isequal(test_fields, {'gamma', 'bound', 'noise'}));
                                        
                    % Get betas for 'full model', gamma-ablated, and bn-ablated
                    betaG = betaExplained{iSub,iPhz}(:,idx_g);
                    betaBN = betaExplained{iSub,iPhz}(:,idx_bn);
                    betaNull = betaExplained{iSub,iPhz}(:,idx_null);
                    betaFull = betaExplained{iSub,iPhz}(:,idx_full);
                    
                    % Inverse of the above, but same semantics. e.g. 'gamma' bar is now beta value
                    % for bn - null ('all except gamma added in' instead of 'gamma removed')
                    [barData(1), errData(1,1), errData(2,1)] = quantileWeight(betaBN-betaNull, 0.68, weights{iSub,iPhz});
                    [barData(2), errData(1,2), errData(2,2)] = quantileWeight(betaG-betaNull, 0.68, weights{iSub,iPhz});
                    [barData(3), errData(1,3), errData(2,3)] = quantileWeight(betaFull, 0.68, weights{iSub,iPhz});
                    [barData(4), errData(1,4), errData(2,4)] = meanci(abb{iSub,iPhz}(:,2), 0.68);
                    
                    h = bar([1 nan], [barData; nan(1,4)]);
                    drawnow; % populate h.xoffset
                    errorbar(h(1).XData(1) + [h.XOffset], barData, errData(1,:), errData(2,:), '.k');
                    h(1).FaceColor = color1;
                    h(2).FaceColor = color2;
                    h(3).FaceColor = color0;
                    h(4).FaceColor = [1 1 1];
                    
                    if iSub*iPhz == 1, legend({'leak (\gamma) = 0', 'bound = \infty', 'full', 'true'}); end
            end
            
            title([subjectIds{iSub} ' :: ' plot_titles{iPhz}]);
            
            ylim([-1 1]);
        end
    end
end

end

function [mu, lo, hi, med] = quantileWeight(data, interval, weight)
weight = weight(:);
mu = sum(data.*weight)/sum(weight);

plo = .5-interval/2;
phi = .5+interval/2;
for j=size(data,2):-1:1
    [vals, isrt] = sort(data(:,j));
    cdf = cumsum(weight(isrt)) / sum(weight);
    lo = interp1q(cdf, vals, plo);
    hi = interp1q(cdf, vals, phi);
    med = interp1q(cdf, vals, 0.5);
end
end