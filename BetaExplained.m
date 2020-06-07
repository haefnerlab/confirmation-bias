function fig_explain = BetaExplained(subjectIds, datadir, memodir)
if nargin < 2, datadir = fullfile(pwd, '..', 'PublishData'); end
if nargin < 3, memodir = fullfile(datadir, '..', 'Precomputed'); end

phase_names = {'hslc', 'lshc', 'both', 'both'};
plot_titles = {'HSLC', 'LSHC', 'Both [HSLC]', 'Both [LSHC]'};
data_idxs = [1 1 1 2];
test_fields = {'gamma', 'bound', 'noise'};

%% Load posterior samples per subject per condition
for iSubject=length(subjectIds):-1:1
    for iPhz=length(phase_names):-1:1
        % Load posterior samples over all parameters
        [samples, fields, ~, ~, params, sigs, choices] = GetITBPosteriorSamples(subjectIds{iSubject}, phase_names{iPhz}, 0, datadir, memodir);
        if strcmp(phase_names{iPhz}, 'both')
            sigs = sigs{data_idxs(iPhz)};
            choices = choices{data_idxs(iPhz)};
            params = params(data_idxs(iPhz));
            
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
        samples = samples(1:50:end, :);
        
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
            
            [meanB, loB, hiB] = meanci(betaExplained{iSubject, iPhz});
            bar(1:length(meanB), meanB, 'FaceColor', 'w');
            errorbar(1:length(meanB), meanB, meanB-loB, hiB-meanB, 'ok');
            
            bar(length(meanB)+1, abb{iSubject,iPhz}(2), 'FaceColor', 'r');
            errorbar(length(meanB)+1, abb{iSubject,iPhz}(2), abb_err{iSubject,iPhz}(2), abb_err{iSubject,iPhz}(2), 'ok');
            
            % Invert naming convention: convert list of ablated fields to list of untouched fields.
            abl = ablations{iPhz};
            allfields = unique(horzcat(abl{:}));
            names = cellfun(@(a) strjoin(setdiff(allfields, a), '+'), abl, 'uniformoutput', false);
            names{cellfun(@isempty, names)} = 'null';
            names = [names {'true'}];
            set(gca, 'XTick', 1:length(names), 'XTickLabel', names);
            xtickangle(60);
            
            title([subjectIds{iSubject} ' :: ' plot_titles{iPhz}]);
            
            ylim([-.3 .3]);
        end
    end
end
end