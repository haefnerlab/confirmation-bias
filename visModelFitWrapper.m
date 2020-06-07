function visModelFitWrapper(truemodel, fitmodel, phase)
% Copied from @ModelComparison
model_info = struct(...
    'name', {'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split', 'ideal'}, ...
    'type', {'is', 'vb-czx', 'itb', 'itb', 'itb', 'itb', 'ideal'}, ...
    'fields', {{'prior_C', 'log_lapse', 'gamma', 'samples'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'step_size', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'neggamma', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma_1', 'gamma_2', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'neggamma_1', 'neggamma_2', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse'}});

model_info = model_info(strcmpi(fitmodel, {model_info.name}));

switch lower(phase)
    case 'lshc'
        phaseid = 2;
        phase = 'LSHC';
    case 'hslc'
        phaseid = 1;
        phase = 'HSLC';
    case 'both'
        phaseid = [2 1];
        phase = 'both';
    otherwise
        error('bad phase');
end

%% Copy logic of ModelComparisonByModel
[params, sigs, choices] = LoadOrRun(@GetGroundTruthSimData, {truemodel, phaseid}, ...
    fullfile('../Precomputed', ['gt-sim-' upper(truemodel) '-' phase '.mat']));
prefix = ['gt-' Model.getModelStringID(params(1), true) '-' lower(phase)];

%% Copy logic of ModelComparison
this_params = params;
fields = model_info.fields;
for iP=1:length(this_params)
    % params may be a struct array
    this_params(iP).model = model_info.type;
end
distribs = Fitting.defaultDistributions(fields, false);
this_params = Fitting.setParamsFields(this_params, fields, cellfun(@(f) distribs.(f).priorrnd(1), fields));
[fits, ~, gpinfo, ~] = LoadOrRun(@Fitting.fitModelBADS, ...
    {this_params, sigs, choices, distribs, struct('prefix', prefix)}, ...
    fullfile('../Precomputed', ['badsfit-' prefix '-' model_info.name '.mat']));

grid_points = gpinfo.x;
loglike = gpinfo.y;
[~, srt] = sort(loglike);
grid_points = grid_points(srt,:);
loglike = loglike(srt);

%% Plot - full model

figure;

ll_vals = cellfun(@(fit) fit.ll, fits);
[~,imax] = max(ll_vals);
for iFit=1:length(fits)
    assert(isequal(fields(:), fits{iFit}(1).fit_fields(:)));
    bestfit(iFit,:) = Fitting.getParamsFields(fits{iFit}, fields);
end

if contains(fitmodel, truemodel, 'ignorecase', true)
    gt = Fitting.getParamsFields(params, fields);
else
    gt = nan(size(bestfit, 2));
end
nF = length(fields);
for iF=1:nF
    for jF=1:nF
        subplot(nF, nF, (iF-1)*nF+jF); hold on;
        if iF==jF
            histogram(grid_points(:,iF), 50);
            yl = ylim;
            for iFit=1:length(fits)
                plot(bestfit(iFit,jF)*[1 1], yl, '-k', 'LineWidth', 2);
            end
            plot(bestfit(imax,jF)*[1 1], yl, '--r', 'LineWidth', 2);
            plot(gt(jF)*[1 1], yl, '-g', 'LineWidth', 2);
        else
            sz = 1+30./(1+exp(-zscore(loglike)));
            c = max(loglike, median(loglike));
            scatter(grid_points(:,jF), grid_points(:,iF), sz, c, 'filled');
            yl = ylim; xl = xlim;
            for iFit=1:length(fits)
                plot(bestfit(iFit,jF)*[1 1], yl, '-k', 'LineWidth', 2, 'HandleVisibility', 'off');
                plot(xl, bestfit(iFit,iF)*[1 1], '-k', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            plot(bestfit(imax,jF)*[1 1], yl, '--r', 'LineWidth', 2, 'DisplayName', 'Best Fit');
            plot(xl, bestfit(imax,iF)*[1 1], '--r', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot(gt(jF)*[1 1], yl, '-g', 'LineWidth', 2, 'DisplayName', 'Ground Truth');
            plot(xl, gt(iF)*[1 1], '-g', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
        if iF==2 && jF == 1, legend(); end
        if iF==nF, xlabel(fields{jF}); end
        if jF==1, ylabel(fields{iF}); end
    end
end

%% Load data - Integrator
% subthresh = true;
% prefix = [subjectId '-' num2str(kernel_kappa) '-' SubjectData.phase '-' num2str(subthresh)];
% [bestfit, train_ll, test_ll, ~, model_names] = IntegratorModelComparison(sigs, choices, 0.98, 100, prefix);

%% Stats - Integrator

% idx_ideal = strcmpi(model_names, 'ideal');
% idx_itb = strcmpi(model_names, 'itb');
% idx_gamma = strcmpi(model_names, 'gamma');
% idx_both = strcmpi(model_names, 'itb_gamma');
%
% train_ll_diff = mean(train_ll - train_ll(idx_ideal, :), 2);
% train_ll_diff_err = std(train_ll - train_ll(idx_ideal, :), [], 2) / sqrt(size(train_ll, 2));
% test_ll_diff = mean(test_ll - test_ll(idx_ideal, :), 2);
% test_ll_diff_err = std(test_ll - test_ll(idx_ideal, :), [], 2) / sqrt(size(test_ll, 2));
%
% [itb_crossings, itb_pk, itb_choices] = integratorStats(bestfit(idx_itb, :), sigs);
% [~, gamma_pk, gamma_choices] = integratorStats(bestfit(idx_gamma, :), sigs);
% [both_crossings, both_pk, both_choices] = integratorStats(bestfit(idx_both, :), sigs);
%
% [true_pk, ~, true_pk_err] = CustomRegression.PsychophysicalKernel(sigs, choices==+1, 1, 0, 500, 1);
%
% pmPlotOptions = struct;
% pmPlotOptions.plotData       = false;
% pmPlotOptions.plotAsymptote  = false;
% pmPlotOptions.plotThresh     = false;
% pmPlotOptions.CIthresh       = false;
% pmPlotOptions.yLabel = 'Percent Chose Left';
% switch phase
%     case 1
%         pmStairVar = 'true_ratio';
%         pmPlotOptions.xLabel = 'True Ratio';
%         [pm_fit, ~, ~, ~] = LoadOrRun(@GaborPsychometric, {SubjectData, 1}, ...
%             fullfile(memodir, ['PM-true_ratio-signed-' SubjectData.subjectID '.mat']));
%     case 2
%         pmPlotOptions.xLabel = 'Noise (\kappa)';
%         [pm_fit, ~, ~, ~] = LoadOrRun(@GaborPsychometric, {SubjectData, -2}, ...
%             fullfile(memodir, ['PM-noise-signed-' SubjectData.subjectID '.mat']));
%         pmStairVar = 'sign_noise';
% end
%
% %% Plots
%
% display_names = cellfun(@(nm) strrep(nm, '_', ' '), model_names, 'uniformoutput', false);
%
% fig = figure;
% subplot(1,5,1); hold on;
% h = bar([train_ll_diff(:) test_ll_diff(:)]'); drawnow;
% errorbar(1+[h.XOffset], train_ll_diff, train_ll_diff_err, 'ok');
% errorbar(2+[h.XOffset], test_ll_diff, test_ll_diff_err, 'ok');
% legend(display_names, 'location', 'best');
% set(gca, 'XTick', [1 2], 'XTickLabel', {'train', 'test'});
% ylabel('\Delta LL per trial');
% title('Model Comparison');
%
%
% subplot(1,5,2); hold on;
% plot(itb_crossings', 'Color', [h(idx_itb).FaceColor .2]);
% plot(both_crossings', 'Color', [h(idx_both).FaceColor .2]);
% errorbar(1:10, mean(itb_crossings, 1), std(itb_crossings, [], 1)/sqrt(size(itb_crossings,1)), 'Color', h(idx_itb).FaceColor);
% errorbar(1:10, mean(both_crossings, 1), std(both_crossings, [], 1)/sqrt(size(both_crossings,1)), 'Color', h(idx_both).FaceColor);
% xlabel('frame');
% ylabel('% bound crossings');
% title('Inferred bound crossings');
%
% subplot(1,5,3); hold on;
% errorbar(true_pk, true_pk_err, '-k', 'LineWidth', 2, 'DisplayName', 'MLE');
% errorbar(mean(itb_pk), std(itb_pk)/sqrt(size(itb_pk,1)), 'Color', h(idx_itb).FaceColor);
% errorbar(mean(gamma_pk), std(gamma_pk)/sqrt(size(gamma_pk,1)), 'Color', h(idx_gamma).FaceColor);
% errorbar(mean(both_pk), std(both_pk)/sqrt(size(both_pk,1)), 'Color', h(idx_both).FaceColor);
% title('PKs');
%
% subplot(1,5,4); hold on;
% plot([bestfit(idx_gamma,:).gamma], [bestfit(idx_itb,:).bound], 'ok');
% plot([bestfit(idx_both,:).gamma], [bestfit(idx_both,:).bound], '+', 'Color', h(idx_both).FaceColor);
% xlabel('gamma');
% ylabel('bound');
% legend({'separate', 'combined'}, 'location', 'best');
% title('Inferred gamma and bound');
%
% [uStim, ~, idxBwd] = unique(SubjectData.(pmStairVar)(trials));
% for iStim=length(uStim):-1:1
%     sim_ch = double(itb_choices(iStim==idxBwd, :));
%     [pm_itb(iStim,1), pm_itb(iStim,2:3)] = binofit(sum(sim_ch(:)), numel(sim_ch), .33);
%     sim_ch = double(gamma_choices(iStim==idxBwd, :));
%     [pm_gamma(iStim,1), pm_gamma(iStim,2:3)] = binofit(sum(sim_ch(:)), numel(sim_ch), .33);
%     sim_ch = double(both_choices(iStim==idxBwd, :));
%     [pm_both(iStim,1), pm_both(iStim,2:3)] = binofit(sum(sim_ch(:)), numel(sim_ch), .33);
% end
%
% subplot(1,5,5); hold on;
% plotPsych(pm_fit, pmPlotOptions);
% errorbar(uStim, pm_itb(:,1), pm_itb(:,1)-pm_itb(:,2), pm_itb(:,3)-pm_itb(:,1), 'Color', h(idx_itb).FaceColor);
% errorbar(uStim, pm_gamma(:,1), pm_gamma(:,1)-pm_gamma(:,2), pm_gamma(:,3)-pm_gamma(:,1), 'Color', h(idx_gamma).FaceColor);
% errorbar(uStim, pm_both(:,1), pm_both(:,1)-pm_both(:,2), pm_both(:,3)-pm_both(:,1), 'Color', h(idx_both).FaceColor);
end