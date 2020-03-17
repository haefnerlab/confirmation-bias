%% Global setup

clear;
rng('shuffle');
true_params = Model.newModelParams('model', 'itb', ...
    'var_x', 0.1, ...
    'gamma', .05, ...
    'save_dir', 'tmp', ...
    'trials', 800, ...
    'updates', 5, ...
    'step_size', 0.01, ...
    'bound', 1, ...
    'noise', .01, ...
    'temperature', .1, ...
    'lapse', 1e-3, ...
    'seed', randi(1e9));
fittable_parameters = {'prior_C', 'lapse', 'gamma', 'sensor_noise', 'var_x', 'noise', 'temperature', 'bound', 'updates', 'samples'};
distribs = Fitting.defaultDistributions(fittable_parameters);

% Simulate ground-truth results
data = Model.genDataWithParams(true_params);
results = Model.runVectorized(true_params, data);

% Model.plotPK(true_params, [1 0 0]);

%% MH sample parameters from their priors and visualize

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
[~, samples, fields] = Fitting.fitChoicesMH(EmptyData, emptyParams, distribs, nSamples, 1, nSamples);

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

%% Investigate effect of # inner-loop iterations on the likelihood

field = 'prior_C';
domain = linspace(0, 1);
islog = false;
prior_info = struct(field, distribs.(field));

inners = [1 10 20];
repeats = 3;

test_params = true_params;

figure;
for i=1:length(inners)
    for j=1:repeats
        log_post = marginallogposterior(@Fitting.choiceModelLogProb, ...
            {prior_info, data, results.choices==+1, inners(i)}, test_params, field, domain);
        post_prob = exp(log_post)/sum(exp(log_post));
        
        subplot(2, length(inners), i); hold on;
        plot(domain, log_post, 'LineWidth', 2);
        if islog, set(gca, 'XScale', 'log'); end
        xlim([min(domain) max(domain)]);
        yl = ylim;
        if j == 1
            plot([test_params.(field) test_params.(field)], yl, '--r');
        end
        title([field ' log post n_{inner}=' num2str(inners(i))]);
        
        subplot(2, length(inners), i+length(inners)); hold on;
        plot(domain, post_prob, 'LineWidth', 2);
        if islog, set(gca, 'XScale', 'log'); end
        xlim([min(domain) max(domain)]);
        if j == 1
            plot([test_params.(field) test_params.(field)], [0 max(post_prob)], '--r');
        end
        title([field ' posteriors n_{inner}=' num2str(inners(i))]);
        drawnow;
    end
end
sgtitle(strrep(Model.getModelStringID(true_params), '_', ' '));

%% MAP inference with BADS

test_params = true_params;
fields = {'prior_C', 'gamma', 'temperature', 'bound', 'lapse'};
logplot = [0 0 1 0 1];
nF = length(fields);
prior_info = struct();
for iF=1:nF
    prior_info.(fields{iF}) = distribs.(fields{iF});
end
x0 = cellfun(@(f) test_params.(f), fields);
nInner = 20;
ll_to_nll = -1;
extra_args = {@Fitting.choiceModelLogProb, test_params, fields, {prior_info, data, results.choices==+1, nInner}, ll_to_nll};

LB = cellfun(@(f) prior_info.(f).lb, fields);
UB = cellfun(@(f) prior_info.(f).ub, fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

opts = bads('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'final';
for iRun=10:-1:1
    x0 = PLB + rand(size(PLB)) .* (PUB - PLB);
    [BESTFIT(iRun,:), ~, EXITFLAG(iRun)] = bads(@loglikefn_wrapper, x0, LB, UB, PLB, PUB, [], opts, extra_args{:});
    for iF=1:nF
        for jF=iF:nF
            subplot(nF, nF, nF*(jF-1)+iF); cla; hold on;
            if iF == jF
                if logplot(iF)
                    histogram(BESTFIT(iRun:end, iF), logspace(log10(PLB(iF)), log10(PUB(iF)), 10));
                    plot(true_params.(fields{iF})*[1 1], [0 4], '--r');
                    set(gca, 'XScale', 'log');
                else
                    histogram(BESTFIT(iRun:end, iF), linspace(PLB(iF), PUB(iF), 10));
                    plot(true_params.(fields{iF})*[1 1], [0 4], '--r');
                end
                title(fields{iF});
            else
                plot(BESTFIT(iRun:end, iF), BESTFIT(iRun:end, jF), 'xk');
                plot(true_params.(fields{iF}), true_params.(fields{jF}), 'or');
                % xlim([PLB(iF), PUB(iF)]);
                % ylim([PLB(jF), PUB(jF)]);
                if logplot(iF), set(gca, 'XScale', 'log'); end
                if logplot(jF), set(gca, 'YScale', 'log'); end
            end
        end
    end
    drawnow;
end

%% Inference with VBMC

test_params = true_params;
fields = {'prior_C', 'gamma', 'temperature', 'bound', 'lapse'};
logplot = [0 0 1 0 1];
nF = length(fields);
prior_info = struct();
for iF=1:nF
    prior_info.(fields{iF}) = distribs.(fields{iF});
end
x0 = cellfun(@(f) test_params.(f), fields);
nInner = 20;
extra_args = {@Fitting.choiceModelLogProb, test_params, fields, {prior_info, data, results.choices==+1, nInner}};

LB = cellfun(@(f) prior_info.(f).lb, fields);
UB = cellfun(@(f) prior_info.(f).ub, fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

opts = vbmc('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'iter';
[VP, ELBO, ELBO_SD, EXITFLAG] = vbmc(@loglikefn_wrapper, x0, LB, UB, PLB, PUB, opts, extra_args{:});

%% VBMC plot

xtrue = cellfun(@(f) test_params.(f), fields);
Xsamp = vbmc_rnd(VP, 1e5);
[fig, ax] = cornerplot(Xsamp, fields, xtrue);

%% Try fitting a subject

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

base_params = true_params;
[params_set, stim_set, choice_set, trial_set] = SubjectDataToModelParams(SubjectData, sigs, kernel_kappa, sensor_noise, base_params);
nonempty = ~cellfun(@isempty, choice_set);
params_set = params_set(nonempty);
stim_set   = stim_set(nonempty);
choice_set = choice_set(nonempty);
trial_set = trial_set(nonempty);

% Visualize transformation of signals from subject data.. set 'signed=true' to view bimodals, or
% signed=false to view just positive lobes.
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
fields = {'prior_C', 'gamma', 'noise', 'bound', 'lapse'};
logplot = [0 1 0 1];
nF = length(fields);
prior_info = struct();
for iF=1:nF
    prior_info.(fields{iF}) = distribs.(fields{iF});
end
nInner = 20;
extra_args = {@Fitting.choiceModelLogProb, params_set, fields, {prior_info, stim_set, choice_set, nInner}};

LB  = cellfun(@(f) prior_info.(f).lb,  fields);
UB  = cellfun(@(f) prior_info.(f).ub,  fields);
PLB = cellfun(@(f) prior_info.(f).plb, fields);
PUB = cellfun(@(f) prior_info.(f).pub, fields);

opts = vbmc('defaults');
opts.UncertaintyHandling = true;
opts.Display = 'final';
for iRun=5:-1:1
    trunstart = tic;
    [VP{iRun}, ELBO(iRun), ELBO_SD(iRun), EXITFLAG(iRun)] = vbmc(@loglikefn_wrapper, [], LB, UB, PLB, PUB, opts, extra_args{:});
    toc(trunstart);
end

[flag, VP_bestfit, best_idx, stats] = vbmc_diagnostics(VP);

%% VBMC plot posteriors
Xsamp = vbmc_rnd(VP{best_idx}, 1e5);
[fig, ax] = cornerplot(Xsamp, fields);

%% VBMC vsualize LPO integration
figure;
bestfit = vbmc_mode(VP{best_idx});
for iSet=1:length(params_set)
    for iF=1:nF
        params_set(iSet).(fields{iF}) = bestfit(iF);
    end
    res = Model.runVectorized(params_set(iSet), stim_set{iSet});
    
    subplotsquare(length(params_set), iSet); hold on;
    % Plot LPO trajectory
    plot(res.lpo', 'Color', [0 0 0 .25]);
    chosepos = choice_set{iSet} == +1;
    % Overlay markers on endpoints according to subject's actual choice
    plot(11*ones(sum(chosepos), 1), res.lpo(chosepos, end), 'og');
    plot(11*ones(sum(~chosepos), 1), res.lpo(~chosepos, end), 'xr');
    % Show bounds as dashed black lines
    plot([1 11], +params_set(iSet).bound*[1 1], '--k');
    plot([1 11], -params_set(iSet).bound*[1 1], '--k');
    % Title according to SI and CI of this group
    title(sprintf('SI=%.2f  CI=%.2f', params_set(iSet).sensory_info, params_set(iSet).category_info));
end

%% Helper functions

function [log_posts] = marginallogposterior(loglikefn, args, params, fieldname, values)
log_posts = arrayfun(@(v) loglikefn(setfield(params, fieldname, v), args{:}), values);
log_posts = log_posts - max(log_posts);
end

function [val] = loglikefn_wrapper(xval, loglikefn, params, fields, args, sgn)
for iPara=1:length(params)
    for i=1:length(fields)
        params(iPara).(fields{i}) = xval(i);
    end
end
if ~exist('sgn', 'var'), sgn = +1; end
val = sgn * loglikefn(params, args{:});
end