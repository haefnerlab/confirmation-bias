function [param_set, stim_set, choice_set] = SubjectDataToModelParams(SubjectData, sensory_noise, base_params)

if ~exist('base_params', 'var'), base_params = Model.newModelParams(); end

% Note: assuming 'true_ratio' and 'noise' (aka kappa) are the only parameters changing

unsigned_ratio = max(SubjectData.true_ratio, 1-SubjectData.true_ratio);
allStimParameters = [SubjectData.noise(:) unsigned_ratio(:)];
unqStim = unique(allStimParameters, 'rows');

for iStim=size(unqStim,1):-1:1
    trials = all(allStimParameters == unqStim(iStim, :), 2);
    [LLOs, sensory_LLOs, empirical_sensory_info] = bpg.signalToLLO(...
        SubjectData.ideal_frame_signals(trials, :), [], unqStim(iStim, 2), sensory_noise);
    
    category_info = unqStim(iStim, 2);
    sensory_info = mean(empirical_sensory_info(:));
    
    % Create model params and add some metadata about experimental params
    this_params = base_params;
    this_params.trials        = size(LLOs, 1);
    this_params.frames        = size(LLOs, 2);
    this_params.category_info = category_info;
    this_params.sensory_info  = sensory_info;
    this_params.p_match       = category_info;
    this_params.var_s         = Model.getEvidenceVariance(sensory_info);
    this_params.data_kappa    = unqStim(iStim, 1);
    this_params.data_ratio    = unqStim(iStim, 2);
    
    param_set(iStim) = this_params;
    % Save model inputs and subject choices -- necessary and sufficient info for fitting. To do so,
    % match LLO from 'sensory' distribution only, since p_match<1 introduces degeneracies in the
    % transformation.
    stim_set{iStim} = Model.lloToEvidence(this_params.var_s, sensory_LLOs);
    choice_set{iStim} = sign(SubjectData.choice(trials) - 0.5);
    
    % DEBUGGING
    colors = lines(size(unqStim,1));
    subplot(1,4,1); hold on;
    sigs = SubjectData.ideal_frame_signals(trials, :);
    plot(sigs(:), stim_set{iStim}(:), '.', 'Color', colors(iStim, :), 'DisplayName', mat2str(unqStim(iStim,:)));
    
    subplot(1,4,2); hold on;
    plot(sigs(:), sensory_LLOs(:), '.', 'Color', colors(iStim, :), 'DisplayName', mat2str(unqStim(iStim,:)));
    
    subplot(1,4,3); hold on;
    plot(sigs(:), empirical_sensory_info(:), '.', 'Color', colors(iStim, :), 'DisplayName', mat2str(unqStim(iStim,:)));
    
    subplot(1,4,4); hold on;
    [f, x] = ksdensity([stim_set{iStim}(:); -stim_set{iStim}(:)]);
    plot(x, f, 'Color', colors(iStim, :));
end

subplot(1,4,1);
xlabel('getSignal()'); ylabel('Model stim');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
legend();

subplot(1,4,2);
xlabel('getSignal()'); ylabel('LLO');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');

subplot(1,4,3);
xlabel('getSignal()'); ylabel('sensory info');

subplot(1,4,4);
xlabel('Model stim');
end