function [fit_params, fit_vals, nll, exitflag] = fitModelToModel(true_params, modeltofit, fields, allow_gamma_neg)
base_params = true_params;
base_params.model = lower(modeltofit);
base_params.allow_gamma_neg = allow_gamma_neg;
distribs = Fitting.defaultDistributions(fields, true, allow_gamma_neg);

data = Model.genDataWithParams(true_params);
true_res = Model.runVectorized(true_params, data);

%% Do fit
LB = cellfun(@(f) distribs.(f).lb, fields);
UB = cellfun(@(f) distribs.(f).ub, fields);
PLB = cellfun(@(f) distribs.(f).plb, fields);
PUB = cellfun(@(f) distribs.(f).pub, fields);

% % DEBUG FIG
% figure;
% subplot(1,2,2); hold on;
% baseline_llo = Model.logLikelihoodOdds(true_params, data);
% baseline_pk = CustomRegression.PsychophysicalKernel(baseline_llo, true_res.choices == +1, 1, 0, 0, 0);
% plot(baseline_pk, '-k', 'LineWidth', 2); yl = ylim;

options = bads('defaults');
options.Display = 'iter';
options.UncertaintyHandling = true;
options.NonlinearScaling = false;

init_vals = rand(size(LB)) .* (PUB - PLB) + PLB;
[fit_vals, nll, exitflag] = bads(...
    @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, x), distribs, data, true_res.choices), ...
    init_vals, LB, UB, PLB, PUB, [], options);

fit_params = Fitting.setParamsFields(base_params, fields, fit_vals);

% % DEBUG FIG
% plot_params = Fitting.setParamsFields(base_params, fields, fit_vals);
% for iF=1:length(fields)
%     for jF=iF:length(fields)
%         subplot(length(fields), 2*length(fields), (jF-1)*2*length(fields)+iF); hold on;
%         scatter(fit_vals(rep, iF), fit_vals(rep, jF), 15, 'filled');
%         % xlim([lb(iF) ub(iF)]);
%         % ylim([lb(jF) ub(jF)]);
%         if iF == 1, ylabel(strrep(fields{jF}, '_', ' ')); end
%         if jF == length(fields), xlabel(strrep(fields{iF}, '_', ' ')); end
%         grid on;
%     end
% end
% res = Model.runVectorized(plot_params, data);
% subplot(1,2,2);
% [pk, ~, pk_err] = CustomRegression.PsychophysicalKernel(baseline_llo, res.choices == +1, 1, 0, 0, 0);
% errorbar(pk, pk_err);
% ylim(yl);
% drawnow;
end