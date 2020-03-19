function [params_samples, samples, fields] = fitChoicesMH(SubjectData, params, distributions, nSamples, maxK, nInner, plotEvery)

fields = fieldnames(distributions);

for iField=1:length(fields)
    if ~isfield(distributions.(fields{iField}), 'logpriorpdf')
        distributions.(fields{iField}).logpriorpdf = ...
            @(x) log(distributions.(fields{iField}).priorpdf(x));
    end
    if ~isfield(distributions.(fields{iField}), 'logproppdf')
        distributions.(fields{iField}).logproppdf = ...
            @(x1, x2) log(distributions.(fields{iField}).proppdf(x1, x2));
    end
end

% TODO - instead of z-scoring, convert experimental data into log-odds
expt_choices = SubjectData.choice;
expt_signals = zeros(size(SubjectData.ideal_frame_signals));
uKappas = unique(SubjectData.noise);
for iKappa=1:length(uKappas)
    trials = SubjectData.noise == uKappas(iKappa);
    sigs = SubjectData.ideal_frame_signals(trials, :);
    expt_signals(trials, :) = sigs / std(sigs(:));
end

    function log_post = logpostpdf(smpl)
        % Fitting.choiceModelLogProb is a stochastic estimator of the log posterior.. average over a few runs
        for i=nInner:-1:1
            log_post(i) = Fitting.choiceModelLogProb(params, distributions, expt_signals, expt_choices, maxK);
        end
        log_post = mean(log_post);
    end

    function x = propose(x)
        for jField=1:length(fields)
            x(jField) = distributions.(fields{jField}).proprnd(x(jField));
        end
    end

    function log_prop = logproppdf(x1, x2)
        log_prop = 0;
        for jField=1:length(fields)
            log_prop = log_prop + distributions.(fields{jField}).logproppdf(x1(jField), x2(jField));
        end
    end

%% Run sampler

initial_values = cellfun(@(f) params.(f), fields);
if nargin < 6
    plotEvery = 100;
end
nBatches = nSamples / plotEvery;
lastSample = initial_values(:)';
samples = [];
for iBatch=1:nBatches
    fprintf('batch %d/%d\t', iBatch, nBatches);
    tic;
    samples(end+1:end+plotEvery,:) = mhsample(lastSample, plotEvery, 'logpdf', ...
        @logpostpdf, 'logproppdf', @logproppdf, 'proprnd', @propose);
    toc;
    lastSample = samples(end, :);
    clf;
    [~, Ax] = plotmatrix(samples);
    for iFild=1:length(fields)
        title(Ax(iField,iField), fields{iField});
    end
    drawnow;
end

%% Assign results back to struct array

for iSample=size(samples, 1):-1:1
    for iField=1:length(fields)
        params_samples(iSample).(fields{iField}) = samples(iSample, iField);
    end
end

end