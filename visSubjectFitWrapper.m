function visSubjectFitWrapper(subjectId, phase)
% model = {'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split', 'ideal'};

kernel_kappa = 0.16;

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

%% Load fit data - full model
prefix = [subjectId '-' num2str(kernel_kappa) '-' num2str(phase)];

for iMod=length(model_info):-1:1
    cache_file = fullfile('../Precomputed', ['qrgfit-' prefix '-' model_info(iMod).name '.mat']);
    if ~exist(cache_file, 'file')
        disp(['Skipping ' cache_file]);
        model_info(iMod).complete = false;
        continue;
    end
    disp(['Loading ' cache_file]);
    model_info(iMod).complete = true;
    ld = load(cache_file);
    model_info(iMod).fit = ld.results{1};
    [~,srt] = sort(ld.results{3}.loglike);
    model_info(iMod).samples = ld.results{2}(srt, :);
    model_info(iMod).loglike = ld.results{3}.loglike(srt);
    if length(phase) == 1
        model_info(iMod).fields = [model_info(iMod).fields {'log_signal_scale'}];
    else
        model_info(iMod).fields = [model_info(iMod).fields {'log_signal_scale_1', 'log_signal_scale_2'}];
    end
end

model_info = model_info([model_info.complete]);

%% Plot - full model
for iMod=1:length(model_info)
    figure;
    mle = Fitting.getParamsFields(model_info(iMod).fit.mle_params, model_info(iMod).fields);
    map = Fitting.getParamsFields(model_info(iMod).fit.map_params, model_info(iMod).fields);
    gp_mle = Fitting.getParamsFields(model_info(iMod).fit.gp_mle_params, model_info(iMod).fields);
    gp_map = Fitting.getParamsFields(model_info(iMod).fit.gp_map_params, model_info(iMod).fields);
    mu = Fitting.getParamsFields(model_info(iMod).fit.mean_params, model_info(iMod).fields);
    nF = length(model_info(iMod).fields);
    for iF=1:nF
        for jF=1:nF
            subplot(nF, nF, (iF-1)*nF+jF); hold on;
            if iF==jF
                histogram(model_info(iMod).samples(:,iF), 50);
                plot(mle(jF)*[1 1], ylim, '-r', 'LineWidth', 2);
                plot(map(jF)*[1 1], ylim, '-y', 'LineWidth', 2);
                plot(gp_mle(jF)*[1 1], ylim, '--r', 'LineWidth', 2);
                plot(gp_map(jF)*[1 1], ylim, '--y', 'LineWidth', 2);
                plot(mu(jF)*[1 1], ylim, '-b', 'LineWidth', 2);
            else
                sz = 1+30./(1+exp(-zscore(model_info(iMod).loglike)));
                c = max(model_info(iMod).loglike, median(model_info(iMod).loglike));
                scatter(model_info(iMod).samples(:,jF), model_info(iMod).samples(:,iF), sz, c, 'filled');                yl = ylim; xl = xlim;
                plot(mle(jF)*[1 1], yl, '-r', 'LineWidth', 2, 'DisplayName', 'MLE');
                plot(xl, mle(iF)*[1 1], '-r', 'LineWidth', 2, 'HandleVisibility', 'off');
                plot(map(jF)*[1 1], yl, '-y', 'LineWidth', 2, 'DisplayName', 'MAP');
                plot(xl, map(iF)*[1 1], '-y', 'LineWidth', 2, 'HandleVisibility', 'off');
                plot(mu(jF)*[1 1], yl, '-b', 'LineWidth', 2, 'DisplayName', 'mean');
                plot(xl, mu(iF)*[1 1], '-b', 'LineWidth', 2, 'HandleVisibility', 'off');
                plot(gp_mle(jF)*[1 1], yl, '--r', 'LineWidth', 2, 'DisplayName', 'GP MLE');
                plot(xl, gp_mle(iF)*[1 1], '--r', 'LineWidth', 2, 'HandleVisibility', 'off');
                plot(gp_map(jF)*[1 1], yl, '--y', 'LineWidth', 2, 'DisplayName', 'GP MAP');
                plot(xl, gp_map(iF)*[1 1], '--y', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            if iF==nF, xlabel(model_info(iMod).fields{jF}); end
            if jF==1, ylabel(model_info(iMod).fields{iF}); end
            if iF==nF && jF==1, legend('location', 'best'); end
        end
    end
    sgtitle([model_info(iMod).name ' :: LL = ' num2str(model_info(iMod).fit.gp_mle_params.ll)]);
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
% errorbar(true_pk, true_pk_err, '-k', 'LineWidth', 2);
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