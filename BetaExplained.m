function fig_explain = BetaExplained(subjectIds, plotstyle, datadir, memodir)
if nargin < 2, plotstyle = 'stacked'; end
if nargin < 3, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 4, memodir = fullfile(datadir, '..', 'Precomputed'); end

% GBN = gamma/bound/noise. inv = inverted.
assert(ismember(lower(plotstyle), {'gbn', 'gbn-inv', 'separate'}));

phase_names = {'hslc', 'lshc', 'both', 'both'};
plot_titles = {'HSLC', 'LSHC', 'Both [HSLC]', 'Both [LSHC]'};
data_idxs = [1 1 1 2];
test_fields = {'gamma', 'bound', 'noise'};

%% Load posterior samples per subject per condition
for iSubject=length(subjectIds):-1:1
    for iPhz=length(phase_names):-1:1
        % Load posterior samples over all parameters
        [samples, fields, ~, ~, params, sigs, choices] = GetITBPosteriorSamples(subjectIds{iSubject}, phase_names{iPhz}, 0, true, datadir, memodir);
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
        thin_idx = round(linspace(1, size(samples,1), min(500, size(samples,1))));
        samples = samples(thin_idx, :);
        
        % Do ablation test on 'test' parameters. Since this involves simulating but not
        % marginalizing over the model simulations, set model to 'itb' rather than 'itb-int' for
        % speed.
        assert(strcmp(params.model, 'itb-int'), 'sanity failed');
        params.model = 'itb';
        [betaExplained{iSubject, iPhz}, ~, ablations{iPhz}] = PKShapeExplained(sigs, params, test_fields, fields, samples);
        
        % Get 'true' slope of exponential PK for this subject x phase
        [abb{iSubject,iPhz}, ~, abb_err{iSubject,iPhz}] = CustomRegression.ExponentialPK(sigs, choices==+1, 1);
    end
end

%% Plot stuff

fig_explain = figure;
for iSubject=1:length(subjectIds)
    for iPhz=1:length(phase_names)
        if ~isempty(betaExplained{iSubject, iPhz})
            subplot(length(subjectIds), length(phase_names), (iSubject-1)*length(phase_names) + iPhz);
            hold on;
            switch lower(plotstyle)
                case 'separate'
                    [meanB, loB, hiB] = meanci(betaExplained{iSubject, iPhz});
                    bar(1:length(meanB), meanB, 'FaceColor', [.9 .9 .9]);
                    errorbar(1:length(meanB), meanB, meanB-loB, hiB-meanB, 'ok');
                    
                    bar(length(meanB)+1, abb{iSubject,iPhz}(2), 'FaceColor', 'r');
                    errorbar(length(meanB)+1, abb{iSubject,iPhz}(2), abb_err{iSubject,iPhz}(2), abb_err{iSubject,iPhz}(2), 'ok');
                    
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
                                        
                    % Null model has all parameters ablated
                    idx_null = cellfun(@(abl) isequal(abl, test_fields), ablations{iPhz});
                    % Full model has all parameters present
                    idx_full = cellfun(@isempty, ablations{iPhz});
                    
                    % Contribution of gamma term quantified by ablating bound and noise
                    idx_g = cellfun(@(abl) isequal(abl, {'bound', 'noise'}), ablations{iPhz});
                    % Contribution of bound+noise quantified by ablating 'gamma'
                    idx_bn = cellfun(@(abl) isequal(abl, {'gamma'}), ablations{iPhz});
                    
                    betaG = betaExplained{iSubject,iPhz}(:,idx_g);
                    betaBN = betaExplained{iSubject,iPhz}(:,idx_bn);
                    betaNull = betaExplained{iSubject,iPhz}(:,idx_null);
                    betaFull = betaExplained{iSubject,iPhz}(:,idx_full);
                    
                    [barData, errData(1,1:4), errData(2,1:4)] = meanci([betaNull betaG betaBN betaFull], 0.95);
                    
                    barData(5) = abb{iSubject, iPhz}(2);
                    errData(:,5) = abb_err{iSubject, iPhz}(2);
                    
                    h = bar([1 nan], [barData; nan(1,5)]);
                    drawnow; % populate h.xoffset
                    errorbar(h(1).XData(1) + [h.XOffset], barData, errData(1,:), errData(2,:), '.k');
                    
                    if iSubject*iPhz == 1, legend({'null', 'gamma', 'bound+noise', 'full', 'true'}); end
                case 'gbn-inv'
                    % In 'inverted' case, contribution of each parameter is the beta that is
                    % lost when that parameter is ablated

                    % This is a targeted analysis comparing gamma/bound/noise... use 'separate' flag
                    % for more generic parameter combinations
                    assert(isequal(test_fields, {'gamma', 'bound', 'noise'}));
                    
                    % Null model has all parameters ablated
                    idx_null = cellfun(@(abl) isequal(abl, test_fields), ablations{iPhz});
                    % Full model has all parameters present
                    idx_full = cellfun(@isempty, ablations{iPhz});
                    
                    % Contribution of gamma term quantified by ablating gamma relative to full
                    idx_g = cellfun(@(abl) isequal(abl, {'gamma'}), ablations{iPhz});
                    % Contribution of bound+noise quantified by ablating b & n relative to full
                    idx_bn = cellfun(@(abl) isequal(abl, {'bound', 'noise'}), ablations{iPhz});
                    
                    betaG = betaExplained{iSubject,iPhz}(:,idx_g);
                    betaBN = betaExplained{iSubject,iPhz}(:,idx_bn);
                    betaNull = betaExplained{iSubject,iPhz}(:,idx_null);
                    betaFull = betaExplained{iSubject,iPhz}(:,idx_full);
                    
                    [barData, errData(1,1:4), errData(2,1:4)] = meanci([betaNull betaFull-betaG betaFull-betaBN betaFull], 0.95);
                    
                    barData(5) = abb{iSubject, iPhz}(2);
                    errData(:,5) = abb_err{iSubject, iPhz}(2);
                    
                    h = bar([1 nan], [barData; nan(1,5)]);
                    drawnow; % populate h.xoffset
                    errorbar(h(1).XData(1) + [h.XOffset], barData, errData(1,:), errData(2,:), '.k');
                    
                    if iSubject*iPhz == 1, legend({'null', 'gamma', 'bound+noise', 'full', 'true'}); end
            end
            
            title([subjectIds{iSubject} ' :: ' plot_titles{iPhz}]);
            
            ylim([-.5 .5]);
        end
    end
end
end