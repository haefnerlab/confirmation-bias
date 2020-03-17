%% Investigate effect of # inner-loop iterations on the likelihood

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

%% Inference with VBMC

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