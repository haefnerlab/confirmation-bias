%% Global setup

clear;
rng('shuffle');
true_params = Model.newModelParams('model', 'itb', ...
    'var_x', 0.1, ...
    'gamma', .05, ...
    'allow_gamma_neg', false, ...
    'save_dir', 'tmp', ...
    'trials', 800, ...
    'updates', 1, ...
    'step_size', 0.01, ...
    'bound', 1, ...
    'noise', .01, ...
    'temperature', .1, ...
    'lapse', 1e-3, ...
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
log_distribs_map = Fitting.defaultDistributions(fittable_parameters, true);
log_distribs_bayes = Fitting.defaultDistributions(log_fittable_parameters, false);

% Simulate ground-truth results
data = Model.genDataWithParams(true_params);
results = Model.runVectorized(true_params, data);

% Model.plotPK(true_params, [1 0 0]);

%% Plot priors over each parameter

figure;
for iPara=1:length(fittable_parameters)
    d = distribs.(fittable_parameters{iPara});
    xs = linspace(d.plb, d.pub, 1001);
    log_ps = d.logpriorpdf(xs);
    ps = exp(log_ps - max(log_ps));
    ps = ps / sum(ps);
    
    checkps = d.priorpdf(xs);
    checkps = checkps / sum(checkps);
    
    subplotsquare(length(fittable_parameters), iPara); hold on;
    plot(xs, ps, '-k', 'LineWidth', 2);
    plot(xs, checkps, '--r', 'LineWidth', 1);
    plot(d.plb*[1 1], ylim, '--b');
    plot(d.pub*[1 1], ylim, '--g');
    title([strrep(fittable_parameters{iPara}, '_', ' ') ' prior']);
    
    if iPara == length(fittable_parameters)
        legend({'exp(logpriorpdf)', 'priorpdf'}, 'location', 'best');
    end
end

clear iPara d xs ps log_ps checkps

%% MH sample parameters from their priors and visualize

error('this block needs updating');

% Sample from the prior by passing empty data
EmptyData = struct('choice', [], 'ideal_frame_signals', [], 'noise', [], 'params', []);
emptyParams = Model.newModelParams('bound', 5, 'lapse', 1e-2, 'gamma', .1, 'var_x', .1, 'noise', .1, 'temperature', .1, 'bound', .5, 'updates', 5, 'samples', 5);
emptyParams = Fitting.sanitize(emptyParams);
emptyParams.sensor_noise = 0.5;

for iF=1:length(fittable_parameters)
    f = fittable_parameters{iF};
    fprintf('%s: begin @ %.1e with logprior %.1e\n', f, emptyParams.(f), ...
        distribs.(f).logpriorpdf(emptyParams.(f)));
end

nSamples = 5000;
[~, samples, fields] = Fitting.fitChoicesMH(EmptyData, emptyParams, distribs, nSamples, 1e3, 1, nSamples);

figure;
nlag = 300;
for i=1:length(fittable_parameters)
    subplot(length(fittable_parameters), 3, 3*(i-1)+[1 2]); hold on;
    plot(samples(:, i));
    yl = ylim;
    plot([1 nSamples], distribs.(fittable_parameters{i}).lb*[1 1], '-r');
    plot([1 nSamples], distribs.(fittable_parameters{i}).ub*[1 1], '-r');
    plot([1 nSamples], distribs.(fittable_parameters{i}).plb*[1 1], '--r');
    plot([1 nSamples], distribs.(fittable_parameters{i}).pub*[1 1], '--r');
    ylim(yl);
    ylabel(fields{i});
    axis tight;
    for j=1:nlag
        lag = j-1;
        acf(j) = corr(samples(1:end-nlag+1, i), samples(1+lag:end-nlag+1+lag, i));
    end
    subplot(length(fittable_parameters), 3, 3*(i-1)+3);
    plot(0:nlag-1, acf);
    ylim([0 1]);
end

%% Investigate effect of max K iterations on the likelihood

field = 'log_temperature';
prior_info = Fitting.defaultDistributions({field}, false);
domain = linspace(prior_info.(field).plb, prior_info.(field).pub);

maxKs = [10 100 200];
repeats = 3;

test_params = true_params;

figure;
for i=1:length(maxKs)
    for j=1:repeats
        log_post(:,j) = arrayfun(...
            @(x) Fitting.choiceModelLogProb(Fitting.setParamsFields(test_params, field, x), prior_info, data, results.choices, maxKs(i)), ...
            domain);
        post_prob(:,j) = log2prob(log_post(:,j));
        
        subplot(2, length(maxKs), i); hold on;
        plot(domain, log_post(:,j), 'LineWidth', 2);
        xlim([min(domain) max(domain)]);
        yl = ylim;
        title([field ' log post K_{max}=' num2str(maxKs(i))]);
        
        subplot(2, length(maxKs), i+length(maxKs)); hold on;
        plot(domain, post_prob(:,j), 'LineWidth', 2);
        xlim([min(domain) max(domain)]);
        title([field ' posteriors K_{max}=' num2str(maxKs(i))]);
        drawnow;
    end
    best_log_post = smooth(mean(log_post, 2), 7, 'rloess');
    subplot(2, length(maxKs), i); hold on;
    plot(domain, best_log_post, '-', 'Color', [0 0 0 .75], 'LineWidth', 2);
    plot(Fitting.getParamsFields(test_params, field)*[1 1], ylim, '--r');
    subplot(2, length(maxKs), i+length(maxKs)); hold on;
    best_post = log2prob(best_log_post);
    plot(domain, best_post/sum(best_post), '-', 'Color', [0 0 0 .75], 'LineWidth', 2);
    plot(Fitting.getParamsFields(test_params, field)*[1 1], ylim, '--r');
    drawnow;
end
sgtitle(strrep(Model.getModelStringID(true_params), '_', ' '));

%% MAP inference fitting model to itself with BADS

test_params = true_params;
fields = {'prior_C', 'gamma', 'log_temperature', 'log_bound', 'log_lapse'};
nF = length(fields);
prior_info = Fitting.defaultDistributions(fields, true, true_params.allow_gamma_neg);
maxK = 200;

LB = cellfun(@(f) prior_info.(f).lb, fields);
UB = cellfun(@(f) prior_info.(f).ub, fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

opts = bads('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'final';
for iRun=10:-1:1
    x0 = PLB + rand(size(PLB)) .* (PUB - PLB);
    [BESTFIT(iRun,:), ~, EXITFLAG(iRun)] = bads(...
        @(x) -Fitting.choiceModelLogProb(Fitting.setParamsFields(test_params, fields, x), prior_info, data, results.choices, maxK), ...
        x0, LB, UB, PLB, PUB, [], opts);
    for iF=1:nF
        for jF=iF:nF
            subplot(nF, nF, nF*(jF-1)+iF); cla; hold on;
            if iF == jF
                histogram(BESTFIT(iRun:end, iF), linspace(PLB(iF), PUB(iF), 10));
                plot(Fitting.getParamsFields(true_params, fields{iF})*[1 1], [0 4], '--r');
                title(strrep(fields{iF}, '_', ' '));
            else
                plot(BESTFIT(iRun:end, iF), BESTFIT(iRun:end, jF), 'xk');
                plot(Fitting.getParamsFields(true_params, fields{iF}), Fitting.getParamsFields(true_params, fields{jF}), 'or');
            end
        end
    end
    drawnow;
end

%% Inference fitting model to itself with VBMC

test_params = true_params;
fields = {'prior_C', 'gamma', 'temperature', 'bound', 'lapse'};
nF = length(fields);
prior_info = struct();
for iF=1:nF
    prior_info.(fields{iF}) = distribs.(fields{iF});
end
x0 = cellfun(@(f) test_params.(f), fields);
maxK = 1e3;
extra_args = {@Fitting.choiceModelLogProb, test_params, fields, {prior_info, data, results.choices, maxK}};

LB = cellfun(@(f) prior_info.(f).lb, fields);
UB = cellfun(@(f) prior_info.(f).ub, fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

opts = vbmc('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'iter';
[VP, ELBO, ELBO_SD, EXITFLAG] = vbmc(@logprobfn_wrapper, x0, LB, UB, PLB, PUB, opts, extra_args{:});

%% Fit model to self VBMC plot

xtrue = cellfun(@(f) test_params.(f), fields);
Xsamp = vbmc_rnd(VP, 1e5);
[fig, ax] = cornerplot(Xsamp, fields, xtrue);

%% Load subject data and convert to model space

% Uppercase constants copied from @PaperFigures
RATIO_PHASE = 1; % aka HSLC
NOISE_PHASE = 2; % aka LSHC
THRESHOLD = 0.7;
DATADIR = fullfile(pwd, '..', 'PublishData');
MEMODIR = fullfile(pwd, '..', 'Precomputed');

suffix = 'noise';
kernel_kappa = 0.16;
sensor_noise = 0; % TODO - for noise > 0 some non-arbitrary decisions still need to be made, such as whether to literally add internal noise, adjust signals' CDF, etc.
subjectId = 'BPGTask-subject07';

SubjectData = LoadAllSubjectData(subjectId, NOISE_PHASE, DATADIR);
sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, kernel_kappa}, ...
    fullfile(MEMODIR, ['perFrameSignals-' subjectId '-' num2str(kernel_kappa) '-noise.mat']));

base_params = Model.newModelParams('model', 'itb', 'gamma', .1, 'noise', .01, 'temperature', .01, 'bound', 3);
[params_set, stim_set, choice_set, trial_set] = SubjectDataToModelParams(SubjectData, sigs, kernel_kappa, sensor_noise, base_params);
nonempty = ~cellfun(@isempty, choice_set);
params_set = params_set(nonempty);
stim_set   = stim_set(nonempty);
choice_set = choice_set(nonempty);
trial_set = trial_set(nonempty);

% Visualize transformation of signals from subject data.. set 'signed=true' to view bimodals, or
% signed=false to view just positive lobes.
figure;
signed = true;
for iSet=1:length(params_set)
    subplotsquare(length(params_set), iSet); cla; hold on;
    this_sigs = sigs(trial_set{iSet}, :);
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

% Prepare model fields, priors, etc for fitting
fields = {'prior_C', 'gamma', 'noise', 'temperature', 'bound', 'lapse'};
nF = length(fields);
empty_prior = struct(); % Use this when just evaluating log likelihood
prior_info = struct(); % Use this to get log posterior
for iF=1:nF
    prior_info.(fields{iF}) = distribs.(fields{iF});
end
maxK = 1e3;

LB  = cellfun(@(f) prior_info.(f).lb,  fields);
UB  = cellfun(@(f) prior_info.(f).ub,  fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

%% Baseline: evaluate subject w.r.t. ideal observer model, fitting only temperature param

ideal_params = arrayfun(@(p) setfield(p, 'model', 'ideal'), params_set);
temps = 0:.1:10;
for iT=length(temps):-1:1
    for iSet=1:length(params_set)
        ideal_params(iSet).temperature = temps(iT);
    end
    ideal_ll(iT) = Fitting.choiceModelLogProb(ideal_params, empty_prior, stim_set, choice_set, maxK);
end

[~,imax] = max(ideal_ll);
% Set to optimal temperature
ideal_params = arrayfun(@(p) setfield(p, 'temperature', temps(imax)), ideal_params);

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
sgtitle({['Ideal Observer behavior on ' subjectId '-translated data'], ['Choice-model LL = ' num2str(ideal_ll(imax))]});

%% Try fitting subject with VBMC

% First, copy best-fit 'temperature' parameter from the ideal observer model.
vbmc_params_set = params_set;
for iSet=1:length(params_set)
    vbmc_params_set(iSet).temperature = ideal_params(iSet).temperature;
end
% Arguments # 2..end to pass to @logprobfn_wrapper
extra_args = {@Fitting.choiceModelLogProb, vbmc_params_set, fields, {prior_info, stim_set, choice_set, maxK}};

opts = vbmc('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'final';
iRun = 1;
flag = -1;
% Fit a minimum of 5 models... then keep going until the vbmc_diagnostics function is satisfied with
% convergence.
while iRun < 5 || flag <= 0
    trunstart = tic;
    [VP{iRun}, ELBO(iRun), ELBO_SD(iRun)] = vbmc(@logprobfn_wrapper, [], LB, UB, PLB, PUB, opts, extra_args{:});
    toc(trunstart);
    if iRun >= 5
        [flag, VP_bestfit, bestRun, stats] = vbmc_diagnostics(VP);
    end
    iRun = iRun+1;
end

%% VBMC :: plot best posteriors
Xsamp = vbmc_rnd(VP{bestRun}, 1e5);
[fig, ax] = cornerplot(Xsamp, fields, [], [LB; UB]);

%% VBMC :: visualize best model's LPO integration
figure;
bestfit = vbmc_mode(VP{bestRun});
for iSet=1:length(vbmc_params_set)
    for iF=1:nF
        vbmc_params_set(iSet).(fields{iF}) = bestfit(iF);
    end
    res = Model.runVectorized(vbmc_params_set(iSet), stim_set{iSet});
    
    subplotsquare(length(vbmc_params_set), iSet); hold on;
    % Plot LPO trajectory
    plot(res.lpo', 'Color', [0 0 0 .25]);
    chosepos = choice_set{iSet} == +1;
    % Overlay markers on endpoints according to subject's actual choice
    plot(11*ones(sum(chosepos), 1), res.lpo(chosepos, end), 'og');
    plot(11*ones(sum(~chosepos), 1), res.lpo(~chosepos, end), 'xr');
    % Show bounds as dashed black lines
    plot([1 11], +vbmc_params_set(iSet).bound*[1 1], '--k');
    plot([1 11], -vbmc_params_set(iSet).bound*[1 1], '--k');
    % Title according to SI and CI of this group
    title(sprintf('SI=%.2f  CI=%.2f', vbmc_params_set(iSet).sensory_info, vbmc_params_set(iSet).category_info));
end
vbmc_ll = Fitting.choiceModelLogProb(params_set, empty_prior, stim_set, choice_set, maxK);
sgtitle({['ITB [VBMC-MAP] behavior on ' subjectId '-translated data'], ['Choice-model LL = ' num2str(vbmc_ll)]});

%% Try fitting subject with BADS

% First, copy best-fit 'temperature' parameter from the ideal observer model.
bads_params_set = params_set;
for iSet=1:length(params_set)
    bads_params_set(iSet).temperature = ideal_params(iSet).temperature;
end
% Arguments # 2..end to pass to @logprobfn_wrapper
extra_args = {@Fitting.choiceModelLogProb, bads_params_set, fields, {prior_info, stim_set, choice_set, maxK}};

opts = bads('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'iter';
% Flip sign of LL --> NLL (bads is a *minimization* tool). Passed as 'sgn' arg to @logprobfn_wrapper
bads_args = extra_args;
bads_args{end+1} = -1;
% Run 10 times with random restarts.. keep the best one.
for iRun=10:-1:1
    trunstart = tic;
    x0 = PLB + rand(size(PLB)).*(PUB-PLB);
    [BESTFIT(iRun,:), NLL(iRun), EXITFLAG(iRun)]  = bads(@logprobfn_wrapper, x0, LB, UB, PLB, PUB, opts, bads_args{:});
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
            title(fields{iF});
            plot(BESTFIT(bestRun,iF)*[1 1], ylim, '--r');
        else
            for iRun=1:length(NLL)
                fcol = 1-(NLL(iRun)-min(NLL)+eps)/(max(NLL)-min(NLL)+eps+1);
                icol = ceil(size(cols,1)*fcol);
                % valid = colored, invalid = white
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
MANUALFIT(:, strcmpi(fields, 'prior_C')) = .5;
MANUALFIT(:, strcmpi(fields, 'gamma')) = 0;
MANUALFIT(:, strcmpi(fields, 'noise')) = 0;
MANUALFIT(:, strcmpi(fields, 'temperature')) = ideal_params.temperature;
MANUALFIT(:, strcmpi(fields, 'bound')) = 40;
MANUALFIT(:, strcmpi(fields, 'lapse')) = 0;
disp(fields);
disp(MANUALFIT);
for iSet=1:length(params_set)
    for iF=1:nF
        bads_params_set(iSet).(fields{iF}) = MANUALFIT(iF);
    end
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
bads_loglike = Fitting.choiceModelLogProb(bads_params_set, empty_prior, stim_set, choice_set, maxK);
bads_logpost = Fitting.choiceModelLogProb(bads_params_set, prior_info, stim_set, choice_set, maxK);
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