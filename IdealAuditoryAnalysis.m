function [] = IdealAuditoryAnalysis(subjectID, directory)
prelimFile = fullfile(directory, 'RawData', [subjectID '-AuditoryDataVolume.mat']);
if ~exist(prelimFile, 'file')
    disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
    disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
    return;
else
    load(prelimFile);
end

order_of_clicks = [squeeze(Preliminary_Data.order_of_clicks(:,1,:)) squeeze(Preliminary_Data.order_of_clicks(:,2,:))];
order_of_clicks = reshape(order_of_clicks, Preliminary_Data.current_trial, []);

%wmle=glmfit(X(:,1:end-1), Y', 'binomial');
subplot(2,4,[1,2,3, 4]); hold on;    % Plot left and right weights for the high contrast/volume case
    colors='br';
    LRWeightsErrors = cell(2,2);
    for LeftRight=[1 2]
        st_idx = (LeftRight-1)*120+1;
        end_idx = st_idx+120-1;
        X = order_of_clicks(:,st_idx:end_idx);
        X=X/std(X(:));
        X=[X ones(size(X,1),1)];  % Add a bias term
        Y = Preliminary_Data.choice(:); % 1 x trials
        wmle=glmfit(X(:,1:end-1), Y, 'binomial');
        %     [wAR1,~,SDebars,~] = autoRegress_logisticAR1_2D(X, ...
        %         Y, [2 Test_Data.number_of_images], 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]',wmle); % Call on Jake's functions
        
        %[wAR1,~,SDebars,~] = autoRegress_logisticRidge(X, ...
        %  Y, [1 Test_Data.number_of_images*2], 0.01, 10.^(0:6)',wmle); % Call on Jake's functions
        
        [wAR1,~,SDebars,~] = autoRegress_logisticAR1(X, Y, 120, 0.01, 10.^(0:6)', [.1 .5 .75 .9 .95 .99]',wmle); % Call on Jake's functions
        wAR1(end)=[];      % Remove bias term
        SDebars(end)=[];
        %SDebars_diff = SDebars;   % Save this for the difference in weights
        LRWeightsErrors{LeftRight,1} = wAR1;
        LRWeightsErrors{LeftRight,2} = SDebars;
        
        % Plot the subject's left weights
        linestyle = ['.-' colors(LeftRight)];
        e = errorbar(wAR1, SDebars, linestyle);
        set(e,'Linewidth',2); hold on;
    end
    legend('Left', 'Right');
%e = errorbar(wAR1(Preliminary_Data.number_of_frames+1:end), SDebars(Preliminary_Data.number_of_frames+1:end), '.-r');  % Red Line
%set(e,'Linewidth',2); hold on;
axis tight
legend('Left weights', 'Right weights','Location','northoutside')

%SDebars = SDebars_diff;
wdiff = LRWeightsErrors{1,1}-LRWeightsErrors{2,1};
%s = SDebars(1:length(wAR1)/2) + SDebars(length(wAR1)/2+1:end);
%e = errorbar([1:Preliminary_Data.number_of_images], w', s', '.-k');
subplot(2,4,[5,6,7,8]); hold on;
e= plot(wdiff, '.-k');
set(e,'Linewidth',2); hold on;
axis tight

Get_Figure('Distribution of clicks in bins');
subplot(2,3,[1,2,3]);
plot((1:Preliminary_Data.number_of_frames),Preliminary_Data.counts(1,:));hold on;
subplot(2,3,[4,5,6]);
plot((1:Preliminary_Data.number_of_frames),Preliminary_Data.counts(2,:));hold on;

end