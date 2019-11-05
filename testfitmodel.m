%% Test MH on PK-defined likelihood (Part I)

% true_params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'trials', 1600, 'updates', 5, 'samples', 5, 'seed', 872810841);
true_params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'trials', 1600, 'updates', 5, 'step_size', 0.01, 'seed', 47429240);
data = Model.genDataWithParams(true_params);
results = Model.runVectorized(true_params, data);

nBoot = 500;
boot_weights = zeros(nBoot, true_params.frames + 1);
parfor iBoot=1:nBoot
    bootIdx = randi([1 true_params.trials], true_params.trials, 1);
    boot_weights(iBoot, :) = CustomRegression.PsychophysicalKernel(data(bootIdx, :), results.choices(bootIdx) == +1, 0, 0, 0);
end
pk_mean = mean(boot_weights, 1);
pk_var = var(boot_weights, [], 1);

true_params.var_s_per_sample = true_params.var_s / true_params.samples;
true_params.alpha = 100;

%% Test MH on PK-defined likelihood (Part II)

distrib = Fitting.defaultDistributions({'var_s_per_sample', 'gamma', 'updates', 'samples', 'lapse'});
samples = Fitting.fitPKsMH(pk_mean, pk_var, init_params, distrib, 1000);

%% Plot samples along with ground truth

fields = fieldnames(samples);
n = length(fields);
for i=1:length(fields)
    ivalues = [samples.(fields{i})];
    truei = true_params.(fields{i});
    
    %% marginal
    subplot(n, n, sub2ind([n n], i, i));
    hold on;
    if all(isinteger(ivalues))
        nbins = length(unique(ivalues));
    else
        nbins = 50;
    end
    [marg, edges] = histcounts(ivalues, nbins, 'Normalization', 'pdf');
    bar((edges(1:end-1)+edges(2:end))/2, marg);
    xvals = linspace(min(ivalues), max(ivalues));
    plot(xvals, arrayfun(distrib.(fields{i}).priorpdf, xvals));
    yl = ylim;
    
    plot([truei truei], yl, '--r');
    
    title(strrep(fields{i}, '_', ' '));
    
    %% pairwise joint
    for j=i+1:length(fields)
        jvalues = [samples.(fields{j})];
        truej = true_params.(fields{j});
        subplot(n, n, sub2ind([n n], i, j));
        hold on;
        plot(ivalues, jvalues, '.');
        plot(truei, truej, 'xr');
        
        if i == 1
            ylabel(strrep(fields{j}, '_', ' '));
        end
        
        if j == n
            xlabel(strrep(fields{i}, '_', ' '));
        end
    end
end

%% Inspect marginal likelihoods of each parameter (PK model)

fields_domains = {'var_s',  logspace(-2, 1), true;
    'var_x',  logspace(-2, 1), true;
    'prior_C',  linspace(0, 1), false;
    'gamma',  linspace(0, 1), false;
    'noise', linspace(0, 2), false;
    'updates',  1:100, false;
    'lapse', linspace(.001, 1), false};

figure;
for i=1:size(fields_domains, 1)
    subplot(1, size(fields_domains, 1), i);
    plotmarginalloglikelihood(@Fitting.pkModelLogLikelihood, {pk_mean, pk_var}, true_params, fields_domains{i, :});
    drawnow;
end

%% Investigate effect of # inner-loop iterations on the likelihood (choice model)

field = 'prior_C';
domain = linspace(0, 1);
islog = false;

inners = [1 10 20];
repeats = 3;

test_params = true_params;
test_params.alpha = 1;

figure;
for i=1:length(inners)
    subplot(1, length(inners), i); hold on;
    for j=1:repeats
        plotmarginalloglikelihood(@Fitting.choiceModelLogLikelihood, ...
            {data, results.choices==+1, inners(i)}, ...
            test_params, field, domain, islog);
        drawnow;
    end
end

%% Inference with VBMC (choice model)

test_params = true_params;
test_params.alpha = 1;
fields = {'prior_C', 'gamma', 'alpha', 'samples'};
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

%% Inspect marginal likelihoods of each parameter (Choice model)

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
    plotmarginalloglikelihood(@Fitting.choiceModelLogLikelihood, {data, results.choices==+1, nInner}, ...
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

subjectId = 'bpgFinaltest-subject09';
SubjectData = LoadAllSubjectData(subjectId, 1);

% TODO - fit scale?
uKappas = unique(SubjectData.noise);
for iKappa=1:length(uKappas)
    trials = SubjectData.noise == uKappas(iKappa);
    sigs = SubjectData.ideal_frame_signals(trials, :);
    SubjectData.ideal_frame_signals(trials, :) = sigs / std(sigs(:));
end

init_params = Model.newModelParams('var_x', .1, 'gamma', .2, 'updates', 2, 'samples', 2, 'var_s', .5, 'p_match', .7);

[~, samples, fields] = Fitting.fitChoicesMH(SubjectData, init_params, distributions, 1000, 100);

[~, Ax] = plotmatrix(samples);
for i=1:length(fields)
    title(Ax(i, i), fields{i});
end

%% Helper functions

function [log_likes] = plotmarginalloglikelihood(loglikefn, args, params, fieldname, values, logx)
log_likes = arrayfun(@(v) loglikefn(setfield(params, fieldname, v), args{:}), values);
log_likes = log_likes - max(log_likes);
plot(values, log_likes, '-k', 'LineWidth', 2);
if logx, set(gca, 'XScale', 'log'); end
xlim([min(values) max(values)]);
yl = ylim;
hold on;
plot([params.(fieldname) params.(fieldname)], yl, '--r');
title(fieldname);
drawnow;
end

function [val] = loglikefn_wrapper(xval, loglikefn, params, fields, args)
for i=1:length(fields)
    params.(fields{i}) = xval(i);
end
val = loglikefn(params, args{:});
end