function [param_set, stim_set, choice_set] = SubjectDataToModelParams(SubjectData, sigs, kernel_kappa, sensory_noise, base_params, si_ci_bin_edges)

if ~exist('base_params', 'var')
    base_params = Model.newModelParams('model', 'is', 'frames', SubjectData.number_of_images, ...
        'gamma', 0.1, 'samples', 5, 'updates', 5, 'var_x', 0.1);
end

% Note: assuming 'true_ratio', 'noise' (aka kappa), and 'contrast' are the only parameters changing

unsigned_ratio = max(SubjectData.true_ratio, 1-SubjectData.true_ratio);
allStimParameters = [SubjectData.noise(:) unsigned_ratio(:) SubjectData.contrast(:)];
unqStim = unique(allStimParameters, 'rows');

% Place an exponential prior on kappa with mean 0.16 (which is approximately the threshold level for
% most subjects)
unqKappas = 0:0.04:0.8;
kappaProbs = exp(-unqKappas / 0.16);
kappa_prior = [unqKappas(:) kappaProbs(:)/sum(kappaProbs)];

stim_gen = struct('kappa', SubjectData.noise(1), 'stim_size', SubjectData.stim_size, 'frames', SubjectData.number_of_images, ...
    'stim_sp_freq_cpp', SubjectData.stim_sp_freq_cpp, 'stim_std_sp_freq_cpp', SubjectData.stim_std_sp_freq_cpp, ...
    'contrast', SubjectData.contrast(1), 'annulus', SubjectData.annulus);

% Debugging only
% figure(1); clf; hold on;
for iStim=size(unqStim,1):-1:1
    trials = all(allStimParameters == unqStim(iStim, :), 2);
    
    sig_sign = sign(SubjectData.frame_categories(trials, :));
    
    stim_gen.contrast = unqStim(iStim, 3);
    [s_hat(trials, :), est_si(trials), kappa_post(trials, :), mu_s_hat] = ...
        bpg.standardizeSignal(sigs(trials, :) .* sig_sign, kernel_kappa, stim_gen, kappa_prior, sensory_noise, false);
    s_hat(trials, :) = s_hat(trials, :) .* sig_sign;

    % Debugging only: figure showing posterior over kappa per unique set of stimuli
    % plot(kappa_prior(:,1), mean(kappa_post(trials, :), 1), 'DisplayName', sprintf('inferred \\kappa (\\kappa=%.2f)', stim_gen.kappa));
    % drawnow;
end

% Debugging only: figure showing inferred sensory information vs raw signal level
% figure(2); hold on;
% plot(mean(sigs .* sign(SubjectData.frame_categories), 2), est_si, '.', 'DisplayName', ['\sigma^2_{int} = ' num2str(sensory_noise)]);
% xlabel({'un-adjusted sig each trial', '(output of bpg.getSignal())'});
% ylabel('effective sensory info');
% legend;

%% Group trials by similar sensory/category info

if ~exist('si_ci_bin_edges', 'var'), si_ci_bin_edges = .5:.05:1; end

bin_centers = (si_ci_bin_edges(1:end-1) + si_ci_bin_edges(2:end))/2;
[ss, cc] = meshgrid(bin_centers);

% Adjust lowest bin so that lo < val <= high inequalities capture all expected data
si_ci_bin_edges(1) = si_ci_bin_edges(1)-eps;

nBins = length(si_ci_bin_edges) - 1;
for sBin=nBins:-1:1
    for cBin=nBins:-1:1
        param_set(cBin, sBin) = base_params;
        param_set(cBin, sBin).sensory_info = ss(cBin, sBin);
        param_set(cBin, sBin).var_s = Model.getEvidenceVariance(ss(cBin, sBin));
        param_set(cBin, sBin).category_info = cc(cBin, sBin);
        param_set(cBin, sBin).p_match = cc(cBin, sBin);
        
        % Find all trials where the ratio was in this 'category info' bin and all inferred/empirical
        % sensory infos are in their corresponding bin.
        ci_trials = SubjectData.ratio > si_ci_bin_edges(cBin) & SubjectData.ratio <= si_ci_bin_edges(cBin+1);
        si_trials = est_si > si_ci_bin_edges(sBin) & est_si <= si_ci_bin_edges(sBin+1);
        tr = ci_trials & si_trials;

        % Record all choices the subject made on these trials
        choice_set{cBin, sBin} = SubjectData.choice(tr);

        % Recording stimuli is less straightforward... whereas 's_hat' all have approximately the
        % same variance per difficulty level with changing means, the model is expecting signals
        % drawn from a gaussian with mean +/- 1 and variance that depends on the sensory info. s_hat
        % is the pseudo-zscored signal, meaning it has variance 1 (adjusted for internal noise) and
        % mean given by mu_s_hat.
        this_s_hat = s_hat(tr, :);
        sign_s_hat = sign(SubjectData.frame_categories(tr))';
        s_hat_mog = mog.create(mu_s_hat, ones(size(mu_s_hat)), mean(kappa_post(tr, :), 1));
        zscore_diff_from_mean = ((this_s_hat .* sign_s_hat) - mog.mean(s_hat_mog)) / sqrt(mog.var(s_hat_mog));
        stim_set{cBin, sBin} = sign_s_hat + sign_s_hat .* zscore_diff_from_mean * sqrt(param_set(cBin, sBin).var_s);
    end
end

nData = cellfun(@length, choice_set);
if sum(nData(:)) ~= length(SubjectData.choice)
    warning('%d of %d trials dropped in conversion process!', ...
        length(SubjectData.choice) - sum(nData(:)), length(SubjectData.choice));
end
end