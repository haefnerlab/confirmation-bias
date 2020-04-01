function [samples, accept, base_params] = sampleIntegratorModel(signals, choices, n_samples, fields, allow_gamma_neg, distribs, varargin)
%% Handle inputs

default_fields = {'prior_C', 'gamma', 'bound', 'log_temperature', 'log_lapse_1', 'log_lapse_2'};
default_values = [.5 0 inf 0 -inf -inf];

if nargin < 4 || isempty(fields), fields = default_fields; end
if nargin < 5 || isempty(allow_gamma_neg), allow_gamma_neg = true; end
if nargin < 6 || isempty(distribs), distribs = Fitting.defaultDistributions(fields, false, allow_gamma_neg); end

unrecognized_fields = setdiff(fields, default_fields);
if ~isempty(unrecognized_fields)
    warning('The following fields are not recognized and will be ignored: {%s}', strjoin(unrecognized_fields, ', '));
end

base_params = struct('allow_gamma_neg', allow_gamma_neg);
base_params = Fitting.setParamsFields(base_params, default_fields, default_values);

%% Set up composite proposal / density functions

% If both 'gamma' and 'log_temperature' are included in this model, their posterior will highly
% correlated since, say, stronger negative gamma leads to higher magnitude integration which need to
% be compensated for by raising the temperature. We therefore adjust the proposal distribution using
% the regressed nearly-linear relationship between them, which has the following slope empirically:
log_temp_per_gamma = -4.8;
prop_adjust_gamma_temp = all(ismember({'gamma', 'log_temperature'}, fields));

    function xnew = proprnd(xold)
        xnew = xold;
        for iField=1:length(fields)
            xnew(iField) = distribs.(fields{iField}).proprnd(xold(iField));
        end
        if prop_adjust_gamma_temp
            % Treat log_temperature proposal as a 'delta' away from the regression line with slope
            % log_temp_per_gamma.
            gIdx = strcmpi(fields, 'gamma');
            ltIdx = strcmpi(fields, 'log_temperature');
            delta_gam = xnew(gIdx) - xold(gIdx);
            xnew(ltIdx) = xnew(ltIdx) + delta_gam * log_temp_per_gamma;
        end
    end

    function logp = logproppdf(xnew, xold)
        logp = 0;
        if prop_adjust_gamma_temp
            % Undo shift of log_temperature in the corresponding block of 'proprnd'
            gIdx = strcmpi(fields, 'gamma');
            ltIdx = strcmpi(fields, 'log_temperature');
            delta_gam = xnew(gIdx) - xold(gIdx);
            xnew(ltIdx) = xnew(ltIdx) - delta_gam * log_temp_per_gamma;
        end
        for iField=1:length(fields)
            if isfield(distribs.(fields{iField}), 'logproppdf')
                logp = logp + distribs.(fields{iField}).logproppdf(xnew(iField), xold(iField));
            else
                logp = logp + log(distribs.(fields{iField}).proppdf(xnew(iField), xold(iField)));
            end
        end
    end

    function lp = log_prior(params)
        lp = 0;
        fields = fieldnames(distribs);
        for iF=1:length(fields)
            lp = lp + distribs.(fields{iF}).logpriorpdf(Fitting.getParamsFields(params, fields{iF}));
        end
    end

    function ll = log_like(params, signals, choices)
        run_results = Model.runIntegratorModel(params, signals);
        ll = sum(log(run_results.prob_choice(choices == +1))) + sum(log(1-run_results.prob_choice(choices ~= +1)));
    end

    function lp = log_post(params, signals, choices)
        lp = log_like(params, signals, choices) + log_prior(params);
    end

%% Do sampling

% Draw 100 points from the prior and initialize to the best one
for iInit=100:-1:1
    init_smpl(iInit, :) = cellfun(@(f) distribs.(f).priorrnd(1), fields);
    init_lp(iInit) = log_post(Fitting.setParamsFields(base_params, fields, init_smpl(iInit,:)), signals, choices);
end
[~,idx] = max(init_lp);
init_smpl = init_smpl(idx, :);

% Run sampling
[samples, accept] = mhsample(init_smpl, n_samples, ...
    'logpdf', @(x) log_post(Fitting.setParamsFields(base_params, fields, x), signals, choices), ...
    'proprnd', @proprnd, ...
    'logproppdf', @logproppdf, ...
    varargin{:});

fprintf('%.1f%% samples accepted\n', 100*accept);

% figure;
% nF = length(fields);
% for i=1:nF
%     for j=i:nF
%         subplot(nF,nF,(j-1)*nF+i); hold on;
%
%         for iRep=1:size(fit_vals, 1)
%             plot(init_vals(iRep, i), init_vals(iRep, j), 'ob');
%             plot(fit_vals(iRep, i), fit_vals(iRep, j), 'xg');
%             plot([init_vals(iRep, i) fit_vals(iRep, i)], [init_vals(iRep, j), fit_vals(iRep, j)], '-k');
%         end
%         xl = xlim; yl = ylim;
%         plot(xl, fit_vals(iMAP, j)*[1 1], '--r');
%         plot(fit_vals(iMAP, i)*[1 1], yl, '--r');
%         xlim(xl); ylim(yl);
%
%         if i==1, ylabel(fields{j}); end
%         if j==nF, xlabel(fields{i}); end
%     end
% end
end