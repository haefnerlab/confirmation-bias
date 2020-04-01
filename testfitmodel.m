%% Global setup

clear;
rng('shuffle');
true_params = Model.newModelParams('model', 'is', ...
    'var_x', 0.1, ...
    'gamma', .05, ...
    'allow_gamma_neg', false, ...
    'save_dir', 'tmp', ...
    'trials', 1000, ...
    'updates', 5, ...
    'samples', 5, ...
    'step_size', 0.01, ...
    'bound', 1, ...
    'noise', .01, ...
    'temperature', .1, ...
    'lapse', 1e-2, ...
    'seed', randi(1e9));

% Get default distributions over all fittable parameters (priors and bounds plotted in the next
% block). 3rd argument to defaultDistributions is whether gamma may be negative.
fittable_parameters = {'prior_C', 'lapse', 'gamma', 'sensor_noise', 'var_x', 'noise', ...
    'temperature', 'bound', 'updates', 'samples'};
distribs = Fitting.defaultDistributions(fittable_parameters, [], false);

% The next lines demonstrate that log(x) can be fit instead of x. The semantics of the prior change
% depending on whether we're doing MAP fitting (e.g. with BADS) or full posterior fitting (e.g. with
% VBMC)
log_fittable_parameters = {'log_lapse', 'log_noise', 'log_temperature', 'log_bound'};
log_distribs_map = Fitting.defaultDistributions(log_fittable_parameters, true);
log_distribs_bayes = Fitting.defaultDistributions(log_fittable_parameters, false);

% Simulate ground-truth results
data = Model.genDataWithParams(true_params);
results = Model.runVectorized(true_params, data);

% Model.plotPK(true_params, [1 0 0]);

%% Plot priors over each parameter

figure;
fields = fieldnames(distribs);
for iPara=1:length(fields)
    d = distribs.(fields{iPara});
    xs = linspace(d.plb, d.pub, 1001);
    log_ps = d.logpriorpdf(xs);
    ps = exp(log_ps - max(log_ps));
    ps = ps / sum(ps);
    
    checkps = d.priorpdf(xs);
    checkps = checkps / sum(checkps);
    
    subplotsquare(length(fields), iPara); hold on;
    plot(xs, ps, '-k', 'LineWidth', 2);
    plot(xs, checkps, '--r', 'LineWidth', 1);
    plot(d.plb*[1 1], ylim, '--b');
    plot(d.pub*[1 1], ylim, '--g');
    title([strrep(fields{iPara}, '_', ' ') ' prior']);
    
    if iPara == length(fields)
        legend({'exp(logpriorpdf)', 'priorpdf'}, 'location', 'best');
    end
end

clear iPara d xs ps log_ps checkps

%% Visualize some (log) posterior marginal slices from the true model

field = 'prior_C';
prior_info = Fitting.defaultDistributions({field}, false);
domain = linspace(prior_info.(field).plb, prior_info.(field).pub, 11);

repeats = 3;

test_params = true_params;

lower_bound = -true_params.trials*log(2);

figure;
for j=1:repeats
    log_post(:,j) = arrayfun(...
        @(x) Fitting.choiceModelLogProb(Fitting.setParamsFields(test_params, field, x), prior_info, data, results.choices), ...
        domain);
    post_prob(:,j) = log2prob(log_post(:,j));

    subplot(1,2,1); hold on;
    plot(domain, log_post(:,j), 'LineWidth', 2);
    xlim([min(domain) max(domain)]);

    subplot(1,2,2); hold on;
    plot(domain, post_prob(:,j), 'LineWidth', 2);
    xlim([min(domain) max(domain)]);
    drawnow;
end

subplot(1,2,1);
best_log_post = smooth(mean(log_post, 2), 7, 'rloess');
plot(domain, best_log_post, '-', 'Color', [0 0 0 .75], 'LineWidth', 2);
plot(Fitting.getParamsFields(test_params, field)*[1 1], ylim, '--r');
plot([min(domain) max(domain)], lower_bound*[1 1], '--k');
title([field ' log posts']);

subplot(1,2,2);
title([field ' posteriors']);
best_post = log2prob(best_log_post);
plot(domain, best_post/sum(best_post), '-', 'Color', [0 0 0 .75], 'LineWidth', 2);
plot(Fitting.getParamsFields(test_params, field)*[1 1], ylim, '--r');

sgtitle(strrep(Model.getModelStringID(true_params), '_', ' '));

%% MAP inference fitting model to itself with BADS

test_params = true_params;
test_params.allow_gamma_neg = true;
fields = {'prior_C', 'gamma', 'log_temperature', 'log_bound', 'log_lapse'};
nF = length(fields);
prior_info = Fitting.defaultDistributions(fields, true, test_params.allow_gamma_neg);

LB = cellfun(@(f) prior_info.(f).lb, fields);
UB = cellfun(@(f) prior_info.(f).ub, fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

opts = bads('defaults');
opts.UncertaintyHandling = true;
opts.NonlinearScaling = false;
opts.Display = 'iter';
for iRun=10:-1:1
    x0 = PLB + rand(size(PLB)) .* (PUB - PLB);
    [BESTFIT(iRun,:), NLL(iRun), EXITFLAG(iRun)] = bads(...
        @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(test_params, fields, x), prior_info, data, results.choices), ...
        x0, LB, UB, PLB, PUB, [], opts);
    for iF=1:nF
        for jF=iF:nF
            subplot(nF, nF, nF*(jF-1)+iF); cla; hold on;
            if iF == jF
                histogram(BESTFIT(iRun:end, iF), linspace(PLB(iF), PUB(iF), 10));
                plot(Fitting.getParamsFields(true_params, fields{iF})*[1 1], [0 4], '--r');
                title(strrep(fields{iF}, '_', ' '));
                xlim([PLB(iF) PUB(iF)]);
            else
                plot(BESTFIT(iRun:end, iF), BESTFIT(iRun:end, jF), 'xk');
                plot(Fitting.getParamsFields(true_params, fields{iF}), Fitting.getParamsFields(true_params, fields{jF}), 'or');
                xlim([PLB(iF) PUB(iF)]);
                ylim([PLB(jF) PUB(jF)]);
            end
        end
    end
    drawnow;
end

%% Load subject data and visualize conversion to model space

% Uppercase constants copied from @PaperFigures
RATIO_PHASE = 1; % aka HSLC
NOISE_PHASE = 2; % aka LSHC
THRESHOLD = 0.7;
DATADIR = fullfile(pwd, '..', 'PublishData');
MEMODIR = fullfile(pwd, '..', 'Precomputed');

subjectId = 'BPGTask-subject07';
SubjectData = LoadAllSubjectData(subjectId, NOISE_PHASE, DATADIR);
kernel_kappa = 0.16;
uid = [subjectId '-' num2str(kernel_kappa) '-' SubjectData.phase];
signals = LoadOrRun(@ComputeFrameSignals, {SubjectData, kernel_kappa}, ...
    fullfile('../Precomputed', ['perFrameSignals-' uid '.mat']));

internal_noise = 1;
[params_set, stim_set, choice_set, trial_set] = SubjectDataToModelParams(SubjectData, signals, kernel_kappa, internal_noise);

% Subselect: low(ish) signals only
lowsigcutoff = 0.85;
islowsig = [params_set.sensory_info] <= lowsigcutoff | [params_set.category_info] <= lowsigcutoff;

% Visualize transformation of signals from subject data.. set 'signed=true' to view bimodals, or
% signed=false to view just positive lobes.
figure;
signed = true;
for iSet=1:length(params_set)
    subplotsquare(length(params_set), iSet); cla; hold on;
    this_sigs = signals(trial_set{iSet}, :);
    if signed
        this_sgn = +1;
    else
        this_sgn = sign(SubjectData.frame_categories(trial_set{iSet}, :));
    end
    plot(this_sigs .* this_sgn, stim_set{iSet} .* this_sgn, '.', 'Color', [0 0 0 .25]);
    xl = xlim;
    yl = ylim;
    [f_old, x_old] = ksdensity(this_sigs(:)  .* this_sgn(:));
    [f_new, x_new] = ksdensity(stim_set{iSet}(:) .* this_sgn(:));
    if signed
        f_new_ideal = mog.pdf(x_new(:), mog.create([-1 +1], sqrt(params_set(iSet).var_s)*[1 1], [.5 .5]));
    else
        f_new_ideal = normpdf(x_new, +1, sqrt(params_set(iSet).var_s));
    end
    plot(x_old, yl(1) + f_old / max(f_old) * diff(yl) / 4, '-r');
    plot(xl(1) + f_new / max(f_new) * diff(xl) / 4, x_new, '-b');
    plot(xl(1) + f_new_ideal / max(f_new_ideal) * diff(xl) / 4, x_new, '-g');
    title(sprintf('SI=%.2f  CI=%.2f', params_set(iSet).sensory_info, params_set(iSet).category_info));
end
sgtitle('Debug Fig: red = sig distrib; blue = transformed distrib; green = target of transformation');

%% Baseline: evaluate subject w.r.t. ideal observer model, fitting only temperature param

ideal_params = arrayfun(@(p) setfield(p, 'model', 'ideal'), params_set);
ideal_fields = {'prior_C', 'log_temperature', 'log_lapse'};
ideal_prior_info = Fitting.defaultDistributions(ideal_fields);

LB  = cellfun(@(f) ideal_prior_info.(f).lb,  ideal_fields);
UB  = cellfun(@(f) ideal_prior_info.(f).ub,  ideal_fields);
PLB = cellfun(@(f) ideal_prior_info.(f).plb, ideal_fields);
PUB = cellfun(@(f) ideal_prior_info.(f).pub, ideal_fields);

% Find best-fit prior, temperature, and lapse rate for the otherwise-ideal model
opts = bads('defaults');
opts.UncertaintyHandling = true;
opts.NonlinearScaling = false;
ideal_bestfit = bads(...
    @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(ideal_params(islowsig), ideal_fields, x), ideal_prior_info, stim_set(islowsig), choice_set(islowsig)), ...
    [0.5 log(5) log(.05)], LB, UB, PLB, PUB, [], opts);
ideal_params = Fitting.setParamsFields(ideal_params, ideal_fields, ideal_bestfit);

emp_p_choice = mean(vertcat(choice_set{islowsig}) == +1);
emp_p_category = mean(sum(vertcat(stim_set{islowsig}),2) > 0);
exp_log_prior_C = log(emp_p_choice/(1-emp_p_choice)) - log(emp_p_category/(1-emp_p_category));
exp_prior_C = 1./(1+exp(-exp_log_prior_C));
fprintf('Expected prior from choice bias = %.3f\n', exp_prior_C);
fprintf('Best fit prior = %.3f\n', ideal_bestfit(1));

%% Evaluate log likelihood of the ideal observer model
[~, ideal_ll_lowsig, var_ideal_ll_lowsig] = Fitting.choiceModelLogProb(ideal_params(islowsig), ideal_prior_info, stim_set(islowsig), choice_set(islowsig));
[~, ideal_ll, var_ideal_ll] = Fitting.choiceModelLogProb(ideal_params, ideal_prior_info, stim_set, choice_set);

% Plot ideal observer behavior on subject data
figure;
for iSet=1:length(ideal_params)
    res = Model.runVectorized(ideal_params(iSet), stim_set{iSet});
    
    subplotsquare(length(ideal_params), iSet); hold on;
    % Plot LPO trajectory scaled by temperature
    plot(res.lpo', 'Color', [0 0 0 .25]);
    chosepos = choice_set{iSet} == +1;
    % Overlay markers on endpoints according to subject's actual choice
    plot(11*ones(sum(chosepos), 1), res.lpo(chosepos, end), 'og');
    plot(11*ones(sum(~chosepos), 1), res.lpo(~chosepos, end), 'xr');
    % Title according to SI and CI of this group
    title(sprintf('SI=%.2f  CI=%.2f', ideal_params(iSet).sensory_info, ideal_params(iSet).category_info));
end
sgtitle({['Ideal Observer behavior on ' subjectId '-translated data'], ...
    sprintf('LL=%.2f +/- %.2f (all)\tLL=%.2f +/- %.2f (low)', ideal_ll, sqrt(var_ideal_ll), ideal_ll_lowsig, sqrt(var_ideal_ll_lowsig))});

%% Plot marginal posterior and marginal likelihood slices for each parameter
figure;
for iPara=1:length(ideal_fields)
    subplot(1, length(ideal_fields), iPara); hold on;
    f = ideal_fields{iPara};
    domain = linspace(ideal_prior_info.(f).lb, ideal_prior_info.(f).ub, 31);
    for iX=length(domain):-1:1
        for iRep=10:-1:1
            [logpost(iX, iRep), loglike(iX, iRep), logvar(iX, iRep), lb] = Fitting.choiceModelLogProb(...
                Fitting.setParamsFields(ideal_params(islowsig), {f}, domain(iX)), ideal_prior_info, stim_set(islowsig), choice_set(islowsig));
        end
    end
    errorbar(domain, mean(logpost, 2), sqrt(mean(logvar, 2))/sqrt(size(logvar,2)));
    errorbar(domain, mean(loglike, 2), sqrt(mean(logvar, 2))/sqrt(size(logvar,2)));
    plot(domain, ideal_prior_info.(f).logpriorpdf(domain), '-k');
    xl = xlim; yl = ylim;
    plot(xl, lb*[1 1], '--k');
    plot(Fitting.getParamsFields(ideal_params(1), f)*[1 1], yl, '--b');
    xlim(xl); ylim(yl);
    legend('log posterior', 'log likelihood', 'log prior', 'lower bound');
    title(strrep(f, '_', ' '));
    drawnow;
end

%% Try fitting ITB to subject with BADS

% Prepare model fields, priors, etc for fitting
fields_to_fit = {'prior_C', 'gamma', 'log_temperature', 'log_bound', 'log_lapse'};
nF = length(fields_to_fit);
prior_info = Fitting.defaultDistributions(fields_to_fit, true, false);

LB  = cellfun(@(f) prior_info.(f).lb,  fields_to_fit);
UB  = cellfun(@(f) prior_info.(f).ub,  fields_to_fit);
PLB = cellfun(@(f) prior_info.(f).plb, fields_to_fit);
PUB = cellfun(@(f) prior_info.(f).pub, fields_to_fit);

opts = bads('defaults');
opts.UncertaintyHandling = true;
opts.NonlinearScaling = false;
opts.Display = 'iter';
% Run 10 times with random restarts.. keep the best one.
for iRun=10:-1:1
    trunstart = tic;
    x0 = PLB + rand(size(PLB)).*(PUB-PLB);
    [BESTFIT(iRun,:), NLL(iRun), EXITFLAG(iRun)]  = bads(...
        @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(params_set(islowsig), fields_to_fit, x), prior_info, stim_set(islowsig), choice_set(islowsig)), ...
        x0, LB, UB, PLB, PUB, [], opts);
    toc(trunstart);
end

%% BADS :: plot parameters across runs

% Select best among models that converged satisfactorily
valid = EXITFLAG > 0;
tmp_ll = -NLL;
tmp_ll(~valid) = -inf;
[~, bestRun] = max(tmp_ll);

figure;
cols = jet;
for iF=1:nF
    for jF=iF:nF
        subplot(nF, nF, nF*(jF-1)+iF); cla; hold on;
        if iF == jF
            histogram(BESTFIT(:, iF), linspace(PLB(iF), PUB(iF), 10));
            title(fields_to_fit{iF});
            plot(BESTFIT(bestRun,iF)*[1 1], ylim, '--r');
        else
            for iRun=1:length(NLL)
                fcol = 1-(NLL(iRun)-min(NLL)+eps)/(max(NLL)-min(NLL)+eps+1);
                icol = ceil(size(cols,1)*fcol);
                % Only best pt is filled in, rest are white
                if iRun == bestRun
                    plot(BESTFIT(iRun, iF), BESTFIT(iRun, jF), 'o', 'Color', cols(icol,:), 'MarkerFaceColor', cols(icol, :));
                else
                    plot(BESTFIT(iRun, iF), BESTFIT(iRun, jF), 'o', 'Color', cols(icol,:), 'MarkerFaceColor', [1 1 1]);
                end
            end
            %             plot(BESTFIT(bestRun,iF), BESTFIT(bestRun,jF), 'xk');
            %             xlim([PLB(iF), PUB(iF)]);
            %             ylim([PLB(jF), PUB(jF)]);
        end
    end
end

%% BADS :: visualize best model's LPO integration
figure;
MANUALFIT = BESTFIT(bestRun, :);
% MANUALFIT(:, strcmpi(fields, 'prior_C')) = .5;
% MANUALFIT(:, strcmpi(fields, 'gamma')) = 0;
% MANUALFIT(:, strcmpi(fields, 'noise')) = 0;
% MANUALFIT(:, strcmpi(fields, 'log_temperature')) = log(ideal_params.temperature);
% MANUALFIT(:, strcmpi(fields, 'log_bound')) = 40;
% MANUALFIT(:, strcmpi(fields, 'log_lapse')) = log(1e-3);
disp(fields_to_fit);
disp(MANUALFIT);
bads_params_set = Fitting.setParamsFields(params_set, fields_to_fit, MANUALFIT);
for iSet=1:length(params_set)
    res = Model.runVectorized(bads_params_set(iSet), stim_set{iSet});
    
    subplotsquare(length(bads_params_set), iSet); hold on;
    % Plot LPO trajectory
    plot(res.lpo', 'Color', [0 0 0 .25]);
    chosepos = choice_set{iSet} == +1;
    % Overlay markers on endpoints according to subject's actual choice
    plot(11*ones(sum(chosepos), 1), res.lpo(chosepos, end), 'og');
    plot(11*ones(sum(~chosepos), 1), res.lpo(~chosepos, end), 'xr');
    yl = ylim;
    % Show bounds as dashed black lines
    plot([1 11], +bads_params_set(iSet).bound*[1 1], '--k');
    plot([1 11], -bads_params_set(iSet).bound*[1 1], '--k');
    ylim(yl);
    % Title according to SI and CI of this group
    title(sprintf('SI=%.2f  CI=%.2f', bads_params_set(iSet).sensory_info, bads_params_set(iSet).category_info));
end
bads_loglike = Fitting.choiceModelLogProb(bads_params_set, empty_prior, stim_set, choice_set);
bads_logpost = Fitting.choiceModelLogProb(bads_params_set, prior_info, stim_set, choice_set);
sgtitle({['ITB behavior on ' subjectId '-translated data'], ['Choice-model LL = ' num2str(bads_loglike) ', LPo = ' num2str(bads_logpost) ', LPri = ' num2str(bads_logpost-bads_loglike)]});

%% Helpers

function p = log2prob(logp, dim)
if nargin < 2
    if size(logp, 1) == 1
        dim = 2;
    else
        dim = 1;
    end
end
p = exp(logp-max(logp, [], dim));
p = p ./ sum(p, dim);
end