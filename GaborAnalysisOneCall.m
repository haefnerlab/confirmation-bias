
function []=GaborAnalysisOneCall(subjectIDs, phase, quit, preliminary, threshold, groupings, directory, ideal_template)
if ~exist('directory','var')
    phase = 0;           % 0 = contrast experiment, 1 = ratio experiment
    directory = pwd;     % Make it equal to the current directory
end

number_of_sessions=length(subjectIDs);
if length(quit) ~=number_of_sessions
    disp('Error in quit and subjectID dimension');
    return;
end
subjectID1=cellstr(subjectIDs{1});
if number_of_sessions>1
    for i=1:length(subjectIDs)
        if quit(i)==1
            SaveQuitDataGabor(subjectIDs{i}, phase, directory);
        end
        if i>1
            subjectID2=cellstr(subjectIDs{i});
            
            subjectID_output=strcat(subjectID1, num2str(i));
            ConcatPrelimGabor(cell2mat(subjectID1), cell2mat(subjectID2), cell2mat(subjectID_output), phase, directory);
            subjectID1=subjectID_output;
        end
    end
end
if preliminary==2
    SaveGaborTestPhaseBelowThreshold(cell2mat(subjectID1), phase, threshold, directory);
end
Gabor_Analysis(cell2mat(subjectID1), groupings, preliminary, phase, directory, ideal_template);

end