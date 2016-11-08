function []=GaborAnalysisOneCall(subjectIDs, phase, quit, preliminary, groupings, directory, ideal_template)
if ~exist('directory','var')
    phase = 0;           % 0 = contrast experiment, 1 = ratio experiment
    directory = pwd;     % Make it equal to the current directory
end
number_of_sessions=length(subjectIDs);
if length(quit) ~=number_of_sessions
    disp('Error in quit and subjectID dimension');
    return;
end
subjectID1=subjectIDs{1};
if number_of_sessions>1
for i=1:length(subjectIDs)
   if quit(i)==1
       SaveQuitDataGabor(subjectIDs{i}, phase, directory);
   end
   if i>1
       subjectID2=subjectIDs{i};
       subjectID_output=strcat(subjectID1, num2str(i));
       data_concat_Gabor(subjectID1, subjectID2, subjectID_output, phase, directory);
       subjectID1=subjectID_output;
   end
end
end
Gabor_Analysis(subjectID1, groupings, preliminary, phase, directory, ideal_template);

end