function [param_set, stim_set, choice_set, trial_set] = SubjectDataToModelParams(SubjectData, sigs, kernel_kappa, sensor_noise, base_params)

if ~exist('base_params', 'var')
    base_params = Model.newModelParams('model', 'is', 'frames', SubjectData.number_of_images, ...
        'gamma', 0.1, 'samples', 5, 'updates', 5, 'var_x', 0.1);
end

% Note: assuming 'true_ratio', 'noise' (aka kappa), and 'contrast' are the only parameters changing

unsigned_ratio = max(SubjectData.true_ratio, 1-SubjectData.true_ratio);
allStimParameters = [SubjectData.noise(:) unsigned_ratio(:) SubjectData.contrast(:)];
[unqStim, ~, idxBwd] = unique(allStimParameters, 'rows');
nUnqStim = size(unqStim, 1);

% Prior on kappa will be a Laplace distribution centered on the true value, as if subjects have some
% but not total information about the true kappa each trial.
unqKappas = (0:0.04:0.8)';
kappaPriorTau = .04;

stim_gen = struct('kappa', SubjectData.noise(1), 'stim_size', SubjectData.stim_size, 'frames', SubjectData.number_of_images, ...
    'stim_sp_freq_cpp', SubjectData.stim_sp_freq_cpp, 'stim_std_sp_freq_cpp', SubjectData.stim_std_sp_freq_cpp, ...
    'contrast', SubjectData.contrast(1), 'annulus', SubjectData.annulus);

% Debugging only
% figure(1); clf; hold on;
% cols = jet(size(unqStim,1));
for iStim=nUnqStim:-1:1
    trials = all(allStimParameters == unqStim(iStim, :), 2);
    
    sig_sign = sign(SubjectData.frame_categories(trials, :));
    
    stim_gen.contrast = unqStim(iStim, 3);
    stim_gen.kappa = unqStim(iStim, 1);
    % Set kappa_prior to peak around the true value
    kappa_prior = [unqKappas, exp(-abs(stim_gen.kappa - unqKappas)/kappaPriorTau)];
    [s_hat(trials, :), est_si(trials), kappa_post(trials, :), mu_s_hat] = ...
        bpg.standardizeSignal(sigs(trials, :) .* sig_sign, kernel_kappa, stim_gen, kappa_prior, sensor_noise, false);
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
% plot(mean(sigs .* sign(SubjectData.frame_categories), 2), est_si, '.', 'DisplayName', ['\sigma^2_{int} = ' num2str(sensor_noise)]);
% xlabel({'un-adjusted sig each trial', '(output of bpg.getSignal())'});
% ylabel('effective sensory info');
% legend;

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
    param_set(iStim).sensory_info = max(.51, min(.99, mean(est_si(tr))));
    param_set(iStim).var_s = Model.getEvidenceVariance(param_set(iStim).sensory_info);
    param_set(iStim).category_info = max(.51, min(.99, unqStim(iStim, 2)));
    param_set(iStim).p_match = param_set(iStim).category_info;
    % Add other Subject's stimulus parameters to the model params for debugging - has no effect on
    % the model behavior
    param_set(iStim).kappa = unqStim(iStim, 1);
    param_set(iStim).contrast = unqStim(iStim, 3);

    % Recording stimuli is less straightforward... whereas 's_hat' all have approximately the
    % same variance per difficulty level with changing means, the model is expecting signals
    % drawn from a gaussian with mean +/- 1 and variance that depends on the sensory info. s_hat
    % is the pseudo-zscored signal, meaning it has variance 1 (adjusted for internal noise) and
    % mean given by mu_s_hat.
    bin_s_hat = s_hat(tr, :);
    bin_sign_s_hat = sign(SubjectData.frame_categories(tr, :));
    bin_kappa_post = kappa_post(tr, :);

    % Each trial, convert signals to model distribution by (i) subtracting the mean, (ii) dividing
    % by standard deviation, (iii) multiplying back in expected standard deviation from the model,
    % and (iv) setting mean to +/- 1 as expected by the model.
    for i=sum(tr):-1:1
        % Each s_hat is assumed to be drawn from a mixture of gaussians (mog), one per kappa,
        % weighted by the inferred posterior over kappa. Thanks to standardization done by
        % bpg.standardizeSignals() above, the 'effective' variance of each mode is 1.
        s_hat_mog = mog.create(mu_s_hat, ones(size(mu_s_hat)), bin_kappa_post(i, :));
        unsigned_centered_s_hat = bin_s_hat(i,:).*bin_sign_s_hat(i,:) - mog.mean(s_hat_mog);
        zscored_s_hat = unsigned_centered_s_hat / sqrt(mog.var(s_hat_mog));
        % Whether the new mean is +/- 1 depends on the *known* category for each frame.
        new_mean = bin_sign_s_hat(i,:);
        stim_set{iStim}(i,:) = new_mean + bin_sign_s_hat(i,:).*zscored_s_hat.*sqrt(param_set(iStim).var_s);
    end
end

nData = cellfun(@length, choice_set);
if sum(nData(:)) ~= length(SubjectData.choice)
    warning('%d of %d trials dropped in conversion process!', ...
        length(SubjectData.choice) - sum(nData(:)), length(SubjectData.choice));
end
end