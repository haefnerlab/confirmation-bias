function [param_set, stim_set, choice_set, trial_set] = SubjectDataToModelParams(SubjectData, kernel_kappa, sensor_noise, base_params, datadir)

if ~exist('datadir', 'var'), datadir=fullfile(pwd, '../RawData'); end
memodir = fullfile(datadir, '../Precomputed');

if ~exist('base_params', 'var') || isempty(base_params)
    base_params = Model.newModelParams('model', 'is', 'frames', SubjectData.number_of_images, ...
        'trials', length(SubjectData.choice), 'gamma', 0.1, 'samples', 5, 'updates', 5, 'var_x', 0.1);
end

% Note: assuming 'true_ratio', 'noise' (aka kappa), and 'contrast' are the only parameters changing

unsigned_ratio = max(SubjectData.true_ratio, 1-SubjectData.true_ratio);
allStimParameters = [SubjectData.noise(:) unsigned_ratio(:) SubjectData.contrast(:)];
[unqStim, ~, idxBwd] = unique(allStimParameters, 'rows');
nUnqStim = size(unqStim, 1);

stim_gen = struct('kappa', SubjectData.noise(1), 'stim_size', SubjectData.stim_size, 'frames', SubjectData.number_of_images, ...
    'stim_sp_freq_cpp', SubjectData.stim_sp_freq_cpp, 'stim_std_sp_freq_cpp', SubjectData.stim_std_sp_freq_cpp, ...
    'contrast', SubjectData.contrast(1), 'annulus', SubjectData.annulus);

sigs = LoadOrRun(@ComputeFrameSignals, {SubjectData, kernel_kappa}, ...
    fullfile(memodir, ['perFrameSignals-' SubjectData.subjectID '-' num2str(kernel_kappa) '-' SubjectData.phase '.mat']));

for iStim=nUnqStim:-1:1
    trials = all(allStimParameters == unqStim(iStim, :), 2);
    
    sig_sign = sign(SubjectData.frame_categories(trials, :));
    
    stim_gen.contrast = unqStim(iStim, 3);
    stim_gen.kappa = unqStim(iStim, 1);
    uid_parts = cellfun(@(f) [f '=' num2str(stim_gen.(f))], fieldnames(stim_gen), 'uniformoutput', false);
    uid = ['bpgstats-' strjoin(uid_parts, '-') '-kernelk=' num2str(kernel_kappa)];
    [mu_sig(iStim), var_sig(iStim), ~, ~] = LoadOrRun(@bpg.getSignalStats, {stim_gen, kernel_kappa}, ...
        fullfile(memodir, [uid '.mat'])); 
end

% Get expected d' of signals at each set of stimulus parameters, accounting for internal noise
dprime = mu_sig ./ sqrt(var_sig + sensor_noise);
[sensory_info, model_var_s] = Model.dprimeToSI(dprime);

%% Construct model params and group trial together per unique stimulus combination
param_set  = repmat(base_params, nUnqStim, 1);
stim_set   = cell(size(param_set));
choice_set = cell(size(param_set));
trial_set  = cell(size(param_set));

for iStim=nUnqStim:-1:1
    tr = idxBwd == iStim;

    param_set(iStim).trials = sum(tr);
    trial_set{iStim} = tr;
    
    % Record all choices the subject made on these trials, using [-1,+1] convention rather than [0,1].
    choice_set{iStim} = sign(double(SubjectData.choice(tr)) - 0.5);
    % By convention, the model expects a column vector of choices
    choice_set{iStim} = choice_set{iStim}(:);

    % Record effective sensory and category information for this set of trials
    param_set(iStim).sensory_info = max(.51, min(.99, sensory_info(iStim)));
    param_set(iStim).var_s = model_var_s(iStim);
    param_set(iStim).category_info = max(.51, min(.99, unqStim(iStim, 2)));
    param_set(iStim).p_match = param_set(iStim).category_info;
    % Add other Subject's stimulus parameters to the model params for debugging - has no effect on
    % the model behavior
    param_set(iStim).kappa = unqStim(iStim, 1);
    param_set(iStim).contrast = unqStim(iStim, 3);

    % Recording stimuli is less straightforward... whereas 's_hat' all have approximately the
    % same variance per difficulty level with changing means, the model is expecting signals
    % drawn from a gaussian with mean +/- 1 and variance that depends on the sensory info.
    bin_sgn = sign(SubjectData.frame_categories(tr, :));
    bin_unsgn_sigs = sigs(tr, :) .* bin_sgn;

    % Z-score signals based on their empirical distribution.
    bin_zscore_sigs = (bin_unsgn_sigs - mu_sig(iStim)) ./ var_sig(iStim);

    % Convert back to magnitudes expected by the model, i.e. Gaussian centered at +1 with variance
    % set by sensory info.
    model_unsgn_sigs = 1 + bin_zscore_sigs.*sqrt(param_set(iStim).var_s);
    stim_set{iStim} = model_unsgn_sigs.*bin_sgn;
end

nData = cellfun(@length, choice_set);
if sum(nData(:)) ~= length(SubjectData.choice)
    warning('%d of %d trials dropped in conversion process!', ...
        length(SubjectData.choice) - sum(nData(:)), length(SubjectData.choice));
end
end