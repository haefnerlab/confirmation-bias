function [params_samples, samples, fields] = fitPKsMH(pk_mean, pk_var, params, distributions, nSamples)
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

    function log_post = logpostpdf(smpl)
        %% prior
        log_prior = 0;
        for jField=1:length(fields)
            % Sum the log priors for each parameter
            log_prior = log_prior + distributions.(fields{jField}).logpriorpdf(smpl(jField));
            
            % Also set each field in 'params' to the current sample value
            params.(fields{jField}) = smpl(jField);
        end
        
        params.var_s = params.var_s_per_sample * params.samples;
        
        %% likelihood
        log_likelihood = Fitting.pkModelLogLikelihood(params, pk_mean, pk_var);
        
        %% Combine to get posterior
        log_post = log_prior + log_likelihood;
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
diagnoseEvery = 100;
nBatches = nSamples / diagnoseEvery;
lastSample = initial_values(:)';
samples = [];
for iBatch=1:nBatches
    fprintf('batch %d/%d\t', iBatch, nBatches);
    tic;
    [samples(end+1:end+diagnoseEvery,:), pAccept] = mhsample(lastSample, diagnoseEvery, 'logpdf', ...
        @logpostpdf, 'logproppdf', @logproppdf, 'proprnd', @propose);
    fprintf('accept %.1f%%\t', 100*pAccept);
    toc;
    lastSample = samples(end, :);
    clf;
    [~, Ax] = plotmatrix(samples);
    for iField=1:length(fields)
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