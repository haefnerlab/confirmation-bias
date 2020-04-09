function [best_params, map_value, ll_train, ll_test] = fitIntegratorModel(train_signals, train_choices, fields, allow_gamma_neg, priors, test_signals, test_choices)
%% Handle inputs

default_fields = {'prior_C', 'gamma', 'bound', 'log_temperature', 'log_lapse1', 'log_lapse2'};
default_values = [.5 0 inf 0 -inf -inf];

if nargin < 3 || isempty(fields), fields = default_fields; end
if nargin < 4 || isempty(allow_gamma_neg), allow_gamma_neg = true; end
if nargin < 5 || isempty(priors), priors = Fitting.defaultDistributions(fields, true, allow_gamma_neg); end

unrecognized_fields = setdiff(fields, default_fields);
if ~isempty(unrecognized_fields)
    warning('The following fields are not recognized and will be ignored: {%s}', strjoin(unrecognized_fields, ', '));
end

base_params = struct('allow_gamma_neg', allow_gamma_neg);
base_params = Fitting.setParamsFields(base_params, default_fields, default_values);

%% Objective function(s)

    function nlp = neg_log_prior(params)
        nlp = 0;
        fields = fieldnames(priors);
        for iF=1:length(fields)
            nlp = nlp - priors.(fields{iF}).logpriorpdf(Fitting.getParamsFields(params, fields{iF}));
        end
    end

    function nll = neg_log_like(params, signals, choices)
        run_results = Model.runIntegratorModel(params, signals);
        nll = -sum(log(run_results.prob_choice(choices == +1))) - sum(log(1-run_results.prob_choice(choices ~= +1)));
    end

    function nlp = neg_log_post(params, signals, choices)
        nlp = neg_log_like(params, signals, choices) + neg_log_prior(params);
    end

%% Fit it

LB = cellfun(@(f) priors.(f).lb, fields);
UB = cellfun(@(f) priors.(f).ub, fields);
PLB = cellfun(@(f) priors.(f).plb, fields);
PUB = cellfun(@(f) priors.(f).pub, fields);

options = bads('defaults');
options.Display = 'iter';
options.NonlinearScaling = false;
for iRep=10:-1:1
    % Sample possible starting points from the prior and start from the best one
    for iInit=10:-1:1
        rand_init(iInit, :) = cellfun(@(f) priors.(f).priorrnd(1), fields);
        nlp(iInit) = neg_log_post(Fitting.setParamsFields(base_params, fields, rand_init(iInit, :)), train_signals, train_choices);
    end
    [~, best_init] = min(nlp);
    init_vals(iRep, :) = rand_init(best_init, :);

    % Fit using BADS
    [fit_vals(iRep, :), nlp(iRep), exitflag(iRep)] = bads(...
        @(x) neg_log_post(Fitting.setParamsFields(base_params, fields, x), train_signals, train_choices), ...
        init_vals(iRep, :), LB, UB, PLB, PUB, [], options);
end

nlp(exitflag <= 0) = inf;
[~, iMAP] = min(nlp);

map_value = -nlp(iMAP);
best_params = Fitting.setParamsFields(base_params, fields, fit_vals(iMAP, :));

if nargout >= 3
    ll_train = -neg_log_like(best_params, train_signals, train_choices) / length(train_choices);
end

if nargout >= 4
    ll_test = -neg_log_like(best_params, test_signals, test_choices) / length(test_choices);
end

% figure;
% nF = length(fields);
% for i=1:nF
%     for j=i:nF
%         subplot(nF,nF,(j-1)*nF+i); hold on;
        
%         for iRep=1:size(fit_vals, 1)
%             plot(init_vals(iRep, i), init_vals(iRep, j), 'ob');
%             plot(fit_vals(iRep, i), fit_vals(iRep, j), 'xg');
%             plot([init_vals(iRep, i) fit_vals(iRep, i)], [init_vals(iRep, j), fit_vals(iRep, j)], '-k');
%         end
%         xl = xlim; yl = ylim;
%         plot(xl, fit_vals(iMAP, j)*[1 1], '--r');
%         plot(fit_vals(iMAP, i)*[1 1], yl, '--r');
%         xlim(xl); ylim(yl);
        
%         if i==1, ylabel(fields{j}); end
%         if j==nF, xlabel(fields{i}); end
%     end
% end
end