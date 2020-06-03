function [aic, mle, ll_err, model_info, fits, sampleses] = ModelComparison(base_params, signals, choices, fit_scale, prefix, model_names, best_so_far, memodir)
if nargin < 7, best_so_far = false; end
if nargin < 8, memodir = fullfile('..', 'Precomputed'); end
% Note: use 'base_params' to set generative model parameters like CI, SI, etc

%% Specify models
% Each struct in the array is a model to fit. Colors are taken from matlab's built-in 'lines', then
% grouped so that 'ideal' is yellow, 'is' and 'vb' are warm colors, and all others are cool colors.
model_info = struct(...
    'name', {'is', 'vb', 'itb', 'itb-gamma', 'itb-split', 'itb-gamma-split', 'ideal'}, ...
    'plotname', {'IS', 'VB', 'ITB^+', 'ITB^±', 'ITB^+-\gamma', 'ITB^±-\gamma', 'ITB^0'}, ...
    'type', {'is', 'vb-czx', 'itb', 'itb', 'itb', 'itb', 'ideal'}, ...
    'fields', {{'prior_C', 'log_lapse', 'gamma', 'samples'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'step_size', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'neggamma', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'gamma_1', 'gamma_2', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse', 'neggamma_1', 'neggamma_2', 'bound', 'log_noise'}, ...
              {'prior_C', 'log_lapse'}}, ...
    'color', {[0.8500 0.3250 0.0980], ...
            [0.6350 0.0780 0.1840], ...
            [0.0000 0.4470 0.7410], ...
            [0.4940 0.1840 0.5560], ...
            [0.4660 0.6740 0.1880], ...
            [0.3010 0.7450 0.9330], ...
            [0.9290 0.6940 0.1250]});

if length(base_params) == 1
    % When fitting only a single condition, remove the 'split' models
    model_info = model_info(cellfun(@(nm) ~contains(nm, 'split'), {model_info.name}));
end

% If supplied, restrict to the ones requested
if nargin >= 5
    names = {model_info.name};
    keep = cellfun(@(nm) any(strcmpi(nm, model_names)), names);
    model_info = model_info(keep);
end

%% Fit each model
use_cache = nargin >= 4 && ~isempty(prefix);
mle = nan(size(model_info));
npara = zeros(size(model_info));
ll_var = nan(size(model_info));
fits = cell(size(model_info));
sampleses = cell(size(model_info));
for iModel=1:length(model_info)
    this_params = base_params;
    if fit_scale
        scale_fields = arrayfun(@(phz) sprintf('log_signal_scale_%d', phz), 1:length(base_params), 'uniformoutput', false);
        fields = [model_info(iModel).fields scale_fields];
    else
        fields = model_info(iModel).fields;
    end
    for iP=1:length(this_params)
        % params may be a struct array
        this_params(iP).model = model_info(iModel).type;
    end
    distribs = Fitting.defaultDistributions(fields, false);
    % Randomly init all params. This is primarily so that Model.isStochastic() returns true when
    % 'noise' is a parameter even if base_params.noise=0. But it also just seems like good practice.
    this_params = Fitting.setParamsFields(this_params, fields, cellfun(@(f) distribs.(f).priorrnd(1), fields));
    if use_cache
        cache_file = fullfile(memodir, ['badsfit-' prefix '-' model_info(iModel).name '.mat']);
        if exist(cache_file, 'file') || ~best_so_far
            [fits{iModel}, sampleses{iModel}, ~, ~] = LoadOrRun(@Fitting.fitModelBADS, ...
                {this_params, signals, choices, distribs, struct('prefix', prefix)}, ...
                cache_file);
        elseif best_so_far
            % Kind of a hack... reconstruct the hash for temporary files used inside @FitModelBADS.
            % Load all the interim search data so far and use the best fit.
            if iscell(signals)
                allsigs = vertcat(signals{:});
                allchoices = vertcat(choices{:});
                input_id = string2hash([this_params(1).model, strjoin(fields), num2str([allsigs(:)' allchoices'])]);
            else
                input_id = string2hash([this_params(1).model, strjoin(fields), num2str([signals(:)' choices'])]);
            end
            eval_checkpoint_files = dir(fullfile('fit-checkpoints', sprintf('%x-bads-eval*.mat', input_id)));
            if ~isempty(eval_checkpoint_files)
                for ieval=length(eval_checkpoint_files):-1:1
                    ld = load(fullfile('fit-checkpoints', eval_checkpoint_files(ieval).name));
                    ld_ll(ieval) = ld.this_optim_params(1).ll;
                    ld_ll_var(ieval) = ld.this_optim_params(1).ll_var;
                end
                [n99,pr_n] = Fitting.bootstrapMinimaRegret(-ld_ll, sqrt(ld_ll_var), 1);
                fprintf('%s [%s] not complete\n\tloaded %d of ~%d interim eval files. Pr[min]=%.1f%%\n', ...
                    prefix, model_info(iModel).name, length(eval_checkpoint_files), n99, 100*pr_n(end));

                % Reload the best one as the fit...
                [~,ibest] = max(ld_ll);
                ld = load(fullfile('fit-checkpoints', eval_checkpoint_files(ibest).name));
                fits{iModel} = {ld.this_optim_params};
                sampleses{iModel} = {};
            else
                % No evals.. panic and return nan
                mle(iModel) = nan;
                ll_var(iModel) = nan;
                fits{iModel} = {};
                sampleses{iModel} = {};
                npara(iModel) = nan;
                continue;
            end
        end
    else
        [fits{iModel}, sampleses{iModel}, ~, ~] = Fitting.fitModelBADS(this_params, signals, choices, distribs, struct('prefix', prefix));
    end
    
    % Select best among all local maxima
    fit_lls = cellfun(@(fit_para) fit_para.ll, fits{iModel});
    [mle(iModel), iBest] = max(fit_lls);
    ll_var(iModel) = fits{iModel}{iBest}.ll_var;
    npara(iModel) = length(fields);
end

%% Compute AIC
valid = ~isnan(mle);
aic(valid) = aicbic(mle(valid), npara(valid));
ll_err(valid) = sqrt(ll_var(valid));
aic(~valid) = nan;
ll_err(~valid) = nan;

end