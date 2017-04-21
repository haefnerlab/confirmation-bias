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

% Add a 'true_ratio' field for analysis.
GaborData.true_ratio = sum(GaborData.frame_categories' > 0, 1) / GaborData.number_of_images;

% Add a 'signed noise' and 'signed contrast' field.
trial_sign = sign(GaborData.true_ratio - 0.5);
GaborData.sign_noise = trial_sign .* GaborData.noise;
GaborData.sign_contrast = trial_sign .* GaborData.contrast;

if ~loaded_one
    error('No data found for %s in %s experiment.', subjectID, expt_type);
end
end