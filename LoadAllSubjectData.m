function [GaborData, sources] = LoadAllSubjectData(subjectID, phase, datadir)
%LOADALLSUBJECTDATA Looks for and concatenates .mat files for this subject
%across all sessions, including Quit sessions.
if nargin < 3, datadir = fullfile(pwd, '..', 'RawData'); end

if phase == 0
    expt_type = 'Contrast';
elseif phase == 1
    expt_type = 'Ratio';
elseif phase == 2
    expt_type = 'Noise';
else
    error('Expected phase 0 for Contrast or 1 for Ratio or 2 for Noise');
end

loaded_one = false;
files = dir(fullfile(datadir, '*.mat'));
sources = {};
for i=1:length(files)
    if startsWith(files(i).name, subjectID)
        if endsWith(files(i).name, ['GaborData' expt_type '.mat'])
            contents = load(fullfile(datadir, files(i).name));
            dataToAppend = contents.GaborData;
        elseif endsWith(files(i).name, ['GaborData' expt_type 'Quit.mat'])
            contents = load(fullfile(datadir, files(i).name));
            if contents.GaborData.current_trial < 10
                continue;
            end
            dataToAppend = TruncateQuitDataGabor(contents.GaborData);
        else
            continue;
        end
        sources = horzcat(sources, files(i).name);
        if ~loaded_one
            GaborData = dataToAppend;
            loaded_one = true;
        else
            GaborData = ConcatGaborData(GaborData, dataToAppend);
        end
    end
end

if ~loaded_one
    error('No data found for %s in %s experiment.', subjectID, expt_type);
end

% Add a 'true_ratio' field for analysis.
GaborData.true_ratio = sum(GaborData.frame_categories' > 0, 1) / GaborData.number_of_images;

% Add a 'signed noise' and 'signed contrast' field.
trial_sign = sign(GaborData.true_ratio - 0.5);
GaborData.sign_noise = trial_sign .* GaborData.noise;
GaborData.sign_contrast = trial_sign .* GaborData.contrast;

% Add a blank 'checksum' field if it wasn't there already
if ~isfield(GaborData, 'checksum')
    GaborData.checksum = zeros(size(GaborData.noise));
end

% Add 'flag_use_old_stimulus_code' field if not already present
if ~isfield(GaborData, 'flag_use_old_stimulus_code')
    GaborData.flag_use_old_stimulus_code = false;
end

%% Compute +/- 1 category means for each signal level and normalize per-frame signals accordingly

% [uKappas, ~, idxKappas] = unique(GaborData.noise);
% for iKappa=1:length(uKappas)
%     kappa = max(uKappas(iKappa), 0.04);
%     savename = sprintf('true_sig_sz%d_sp%.3f_spstd%.3f_kap%.3f_ann%.3f.mat', GaborData.stim_size, ...
%         GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, kappa, GaborData.annulus);
%     [~,~,sigs] = LoadOrRun(@bpg.getTrueSignal, {1000, GaborData.stim_size, GaborData.stim_sp_freq_cpp, ...
%         GaborData.stim_std_sp_freq_cpp, 0, kappa, GaborData.annulus, kappa}, ...
%         fullfile('..', 'Precomputed', savename));
%     scale = mean(sigs);
%     trials = idxKappas == iKappa;
%     GaborData.norm_frame_signals(trials, :) = GaborData.ideal_frame_signals(trials, :) / scale;
% end

end