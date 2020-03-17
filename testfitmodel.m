%% Global setup

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
    'seed', randi(1e9));
fittable_parameters = {'prior_C', 'lapse', 'gamma', 'sensor_noise', 'var_x', 'noise', 'temperature', 'bound', 'updates', 'samples'};
distribs = Fitting.defaultDistributions(fittable_parameters);

% Simulate ground-truth results
data = Model.genDataWithParams(true_params);
results = Model.runVectorized(true_params, data);

Model.plotPK(true_params, [1 0 0]);

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

inners = [1 10 20];
repeats = 3;

test_params = true_params;

figure;
for i=1:length(inners)
    for j=1:repeats
        log_like = marginalloglikelihood(@Fitting.choiceModelLogLikelihood, ...
            {data, results.choices==+1, inners(i)}, test_params, field, domain);
        log_pri = distribs.(field).logpriorpdf(domain);
        post_prob = exp(log_like + log_pri)/sum(exp(log_like + log_pri));
        
        subplot(2, length(inners), i); hold on;
        plot(domain, log_like, 'LineWidth', 2);
        if islog, set(gca, 'XScale', 'log'); end
        xlim([min(domain) max(domain)]);
        yl = ylim;
        if j == 1
            plot([test_params.(field) test_params.(field)], yl, '--r');
        end
        title([field ' LL n_{inner}=' num2str(inners(i))]);
        
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
test_params.temperature = 1;
fields = {'prior_C', 'gamma', 'temperature', 'samples', 'lapse'};
x0 = cellfun(@(f) test_params.(f), fields);

extra_args = {@Fitting.choiceModelLogLikelihood, test_params, fields, {data, results.choices==+1, inners(i)}};

LB = [0 0 0 1 0];
UB = [1 1 100 1000 .5];
PLB = [.3 0 .5 1 1e-6];
PUB = [.7 .5 10 100 .3];

[BESTFIT, ~, EXITFLAG] = bads(@loglikefn_wrapper, [], LB, UB, PLB, PUB, [], extra_args{:});

%% Inference with VBMC

test_params = true_params;
test_params.temperature = 1;
fields = {'prior_C', 'gamma', 'temperature', 'samples'};
x0 = cellfun(@(f) test_params.(f), fields);

extra_args = {@Fitting.choiceModelLogLikelihood, test_params, fields, {data, results.choices==+1, inners(i)}};

LB = [0 0 0 1];
UB = [1 1 100 100];
PLB = LB;
PUB = UB;

vbmc_options = vbmc('defaults');
vbmc_options.UncertaintyHandling = 'yes';
[VP, ELBO, ELBO_SD, EXITFLAG] = vbmc(@loglikefn_wrapper, x0, LB, UB, PLB, PUB, vbmc_options, extra_args{:});

%% VBMC plot

Xsamp = vbmc_rnd(VP, 1e5);
[fig, ax] = cornerplot(Xsamp, fields);

%% Inspect marginal likelihoods of each parameter

% fields_domains = {'var_s_per_sample',  logspace(-2, 1), true;
%     'var_x',  logspace(-2, 1), true;
%     'prior_C',  linspace(0, 1), false;
%     'gamma',  linspace(0, 1), false;
%     'updates',  1:100, false;
%     'samples',  1:100, false;
%     'lapse', linspace(.001, 1), false;

fields_domains = {'prior_C',  linspace(0, 1), false};

nInner = 100;

figure;
for i=1:size(fields_domains, 1)
    subplotsquare(size(fields_domains, 1), i);
    marginalloglikelihood(@Fitting.choiceModelLogLikelihood, {data, results.choices==+1, nInner}, ...
        true_params, fields_domains{i, :});
end

%% Try fitting the model to itself on choices

distributions = Fitting.defaultDistributions();
data = Model.genDataWithParams(true_params);
results = Model.runVectorized(true_params, data);

DummySubjectData = struct(...
    'current_trial', true_params.trials, ...
    'ideal_frame_signals', data, ...
    'choice', results.choices, ...
    'noise', true_params.var_s * ones(true_params.trials, 1), ...
    'ratio', true_params.p_match * ones(true_params.trials, 1));

[samples, num_samples, fields] = Fitting.fitChoicesMH(DummySubjectData, true_params, distributions, 1000, 10);

[~, Ax] = plotmatrix(samples);
for i=1:length(fields)
    title(Ax(i, i), fields{i});
end

%% Try fitting a subject

kernel_kappa = 0.16;
subjectId = 'BPGTask-subject07';
SubjectData = LoadAllSubjectData(subjectId, NOISE_PHASE, DATADIR);
sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, kernel_kappa}, ...
    fullfile(MEMODIR, ['perFrameSignals-' subjectId '-' num2str(kernel_kappa) '-noise.mat']));

base_model = Model.newModelParams('model', 'itb');
[params_set, stim_set, choice_set] = SubjectDataToModelParams(SubjectData, sigs, kernel_kappa, 1, base_model);
nonempty = ~cellfun(@isempty, choice_set);
params_set = params_set(nonempty);
stim_set   = stim_set(nonempty);
choice_set = choice_set(nonempty);

fields = {'bound', 'gamma', 'temperature', 'log_noise', 'log_lapse'};
lb = [0 0 0 -12 -12];
ub = [10 1 10 2 log10(.5)];
plb = [1 0 .3 -6 -3];
pub = [5 .5 1 -1 -1];

init_vals = rand(size(plb)).*(pub-plb) + plb;
bestfit = bads(@(x) -loglikefn_wrapper(x, @Fitting.choiceModelLogLikelihood, params_set, fields, {stim_set, choice_set, 10}), ...
    init_vals, lb, ub, plb, pub);

% Plot PK for fit model (low sig only)
eval_params = params_set(1);
eval_params.trials = 10000;
for i=1:length(fields)
    if startsWith(fields{i}, 'log_')
        eval_params.(fields{i}(5:end)) = 10.^bestfit(i);
    else
        eval_params.(fields{i}) = bestfit(i);
    end
end
Model.plotPK(eval_params, [1 0 100]);

% Plot PK for subject
[pk, ~, pk_err] = CustomRegression.PsychophysicalKernel(vertcat(stim_set{1:4}), horzcat(choice_set{1:4}) == +1, 1, 0, 100, 1);
hold on;
errorbar(pk, pk_err);

%% Helper functions

function [log_likes] = marginalloglikelihood(loglikefn, args, params, fieldname, values)
log_likes = arrayfun(@(v) loglikefn(setfield(params, fieldname, v), args{:}), values);
log_likes = log_likes - max(log_likes);
end

function [val] = loglikefn_wrapper(xval, loglikefn, params, fields, args)
for iPara=1:length(params)
    for i=1:length(fields)
        if startsWith(fields{i}, 'log_')
            params(iPara).(fields{i}(5:end)) = 10.^xval(i);
        else
            params(iPara).(fields{i}) = xval(i);
        end
    end
end
val = loglikefn(params, args{:});
end