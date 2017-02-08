function [GaborData, image_collection, sources] = LoadAllSubjectData(subjectID, phase, datadir)
%LOADALLSUBJECTDATA Looks for and concatenates .mat files for this subject
%across all sessions, including Quit sessions.
if nargin < 3, datadir = fullfile(pwd, '..', 'RawData'); end

if phase == 0
    expt_type = 'Contrast';
elseif phase == 1
    expt_type = 'Ratio';
else
    error('Expected phase 0 for Contrast or 1 for Ratio');
end

loaded_one = false;
files = dir(fullfile(datadir, '*.mat'));
sources = {};
for i=1:length(files)
    if startsWith(files(i).name, subjectID)
        if endsWith(files(i).name, ['GaborData' expt_type '.mat'])
            contents = load(fullfile(datadir, files(i).name));
            dataToAppend = contents.GaborData;
            imagesToAppend = contents.image_collection;
        elseif endsWith(files(i).name, ['GaborData' expt_type 'Quit.mat'])
            contents = load(fullfile(datadir, files(i).name));
            [dataToAppend, imagesToAppend] = TruncateQuitDataGabor(contents.GaborData, contents.image_collection);
        else
            continue;
        end
        sources = horzcat(sources, files(i).name);
        if ~loaded_one
            GaborData = dataToAppend;
            image_collection = imagesToAppend;
            loaded_one = true;
        else
            [GaborData, image_collection] = ConcatGaborData(...
                GaborData, image_collection, ...
                dataToAppend, imagesToAppend);
        end
    end
end

if ~loaded_one
    error('No data found for %s in %s experiment.', subjectID, expt_type);
end
end