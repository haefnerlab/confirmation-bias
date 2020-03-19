function [param_set, stim_set, choice_set, trial_set] = SubjectDataToModelParams(SubjectData, sigs, kernel_kappa, sensory_noise, base_params, si_ci_bin_edges)

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
% cols = jet(size(unqStim,1));
for iStim=size(unqStim,1):-1:1
    trials = all(allStimParameters == unqStim(iStim, :), 2);
    
    sig_sign = sign(SubjectData.frame_categories(trials, :));
    
    stim_gen.contrast = unqStim(iStim, 3);
    stim_gen.kappa = unqStim(iStim, 1);
    % % Set kappa_prior to peak around the true value
    % kappa_prior(:,2) = exp(-abs(stim_gen.kappa - kappa_prior(:,1))/.04);
    [s_hat(trials, :), est_si(trials), kappa_post(trials, :), mu_s_hat] = ...
        bpg.standardizeSignal(sigs(trials, :) .* sig_sign, kernel_kappa, stim_gen, kappa_prior, sensory_noise, false);
    s_hat(trials, :) = s_hat(trials, :) .* sig_sign;

    % Debugging  figure option 1: figure showing true and predicted CDF of s_hat for this set of stimuli
    % tmp_s = s_hat(trials, :) .* sig_sign;
    % tmp_s = sort(tmp_s(:));
    % [~,svals] = ksdensity(tmp_s(:));
    % plot(svals, normcdf(svals, mu_s_hat(abs(unqKappas - unqStim(iStim,1)) < eps), 1), 'Color', cols(iStim,:));
    % plot(tmp_s, linspace(0, 1, length(tmp_s)), 'Color', cols(iStim,:), 'Marker', '.');
    % plot(mu_s_hat(abs(unqKappas - unqStim(iStim,1)) < eps)*[1 1], [0 1], 'Color', cols(iStim,:));
    % this_s_hat = s_hat(trials, :);
    % fprintf('%s\t%s\n', mat2str(unqStim(iStim, :)), num2str(var(this_s_hat(:) .* sig_sign(:))));
    
    % Debugging figure option 2: figure showing posterior over kappa per unique set of stimuli
    % plot(kappa_prior(:,1), mean(kappa_post(trials, :), 1), 'DisplayName', sprintf('inferred \\kappa (\\kappa=%.2f)', stim_gen.kappa));
    % drawnow;
end
% xlim([-4 14]);

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

% Adjust edge bins so that lo < val <= high inequalities capture all expected data
si_ci_bin_edges(1) = si_ci_bin_edges(1)-eps;

% Round estimated sensory info to the nearest 1/1000th to deal with some numerical stability issues
est_si = round(est_si * 1000) / 1000;

nBins = length(si_ci_bin_edges) - 1;
stim_set = cell(nBins, nBins);
choice_set = cell(nBins, nBins);
trial_set = cell(nBins, nBins);

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
        
        param_set(cBin, sBin).trials = sum(tr);
        trial_set{cBin, sBin} = tr;
        
        if ~any(tr), continue; end

        % Record all choices the subject made on these trials
        choice_set{cBin, sBin} = SubjectData.choice(tr);

        % Recording stimuli is less straightforward... whereas 's_hat' all have approximately the
        % same variance per difficulty level with changing means, the model is expecting signals
        % drawn from a gaussian with mean +/- 1 and variance that depends on the sensory info. s_hat
        % is the pseudo-zscored signal, meaning it has variance 1 (adjusted for internal noise) and
        % mean given by mu_s_hat.
        bin_s_hat = s_hat(tr, :);
        bin_sign_s_hat = sign(SubjectData.frame_categories(tr, :));
        bin_kappa_post = kappa_post(tr, :);

        % Each trial, convert signals to model distribution by 'pseudo z-scoring', i.e. (i)
        % subtracting the mean, (ii) dividing by standard deviation, (iii) multiplying back in
        % expected standard deviation from the model, and (iv) setting mean to +/- 1 as expected by
        % the model.
        for i=sum(tr):-1:1
            % Each s_hat is assumed to be drawn from a mixture of gaussians (mog), one per kappa,
            % weighted by the inferred posterior over kappa. Thanks to standardization done by
            % bpg.standardizeSignals() above, the 'effective' variance of each mode is 1.
            s_hat_mog = mog.create(mu_s_hat, ones(size(mu_s_hat)), bin_kappa_post(i, :));
            unsigned_centered_s_hat = bin_s_hat(i,:).*bin_sign_s_hat(i,:) - mog.mean(s_hat_mog);
            zscored_s_hat = unsigned_centered_s_hat / sqrt(mog.var(s_hat_mog));
            % Whether the new mean is +/- 1 depends on the *known* category for each frame.
            new_mean = bin_sign_s_hat(i,:);
            stim_set{cBin, sBin}(i,:) = new_mean + bin_sign_s_hat(i,:).*zscored_s_hat.*sqrt(param_set(cBin, sBin).var_s);
        end
    end
end

nData = cellfun(@length, choice_set);
if sum(nData(:)) ~= length(SubjectData.choice)
    warning('%d of %d trials dropped in conversion process!', ...
        length(SubjectData.choice) - sum(nData(:)), length(SubjectData.choice));
end
end