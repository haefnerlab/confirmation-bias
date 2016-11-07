function [] = IdealAnalysis(subjectID,directory)




prelimFile = [directory 'RawData/' subjectID '-AuditoyDataVolume.mat'];
if ~exist(prelimFile, 'file')
    disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
    disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
    return;
else
    load(prelimFile); % Load Preliminary_Data
end

order_of_clicks = [squeeze(Preliminary_Data.order_of_clicks(:,1,:)) squeeze(Preliminary_Data.order_of_clicks(:,2,:))];
order_of_clicks = reshape(order_of_clicks, Preliminary_Data.current_trial, []);
X = order_of_clicks(:,:);
X=X/std(X(:));
Y = Preliminary_Data.choice(:);
[weights,~ , errors,a,b,c] = CustomRegression.PsychophysicalKernel(X, Y,[1], [0],[100000],true);
f= errorbar(weights(1:Preliminary_Data.number_of_frames), errors(1:Preliminary_Data.number_of_frames),'.-b');
hold on;
e = errorbar(weights(Preliminary_Data.number_of_frames+1:end-1), errors(Preliminary_Data.number_of_frames+1:end-1), '.-r');
axis tight;

end