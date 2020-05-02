clear all
clc
close all
%%
% load synthetic data information
load('../dscData/syntheticDataCB/genDataInfo.mat')
genDataInfo.priorC_list = genDataInfo.prior_C_list;
cond = 3;
resFolder = ['../dscData/resultsSyntheticDataCB'];
resNames = dir([resFolder,'/lshc_res*.mat']);
%%
nRound = 2;
[priorCError,gammaError,lapseError,samplesError] = deal(zeros(length(resNames),1));
for i = 1:length(resNames)
    name = [resFolder,'/',resNames(i).name];
    fieldstofit = {'priorC','gamma','lapse','samples'};
    % extract the ground truth of this dataset
    for k = 1:numel(fieldstofit)
        ind =  strfind(name,fieldstofit{k}) + numel(fieldstofit{k});
        eval(sprintf('paraList = genDataInfo.%s_list;',fieldstofit{k}));
        eval(sprintf('groundTruth.%s = paraList(str2double(name(ind)));',fieldstofit{k}));
    end
    load(name);
    
    
    subplot(2,4,[1:4])
    plot(log{1}.loss);
    subplot(2,4,5);histogram(squeeze(posterSamples(nRound,:,1)));hold on;line(groundTruth.priorC*[1,1],[0,50],'color','red');hold off
    subplot(2,4,6);histogram(squeeze(posterSamples(nRound,:,2)));hold on;line(groundTruth.gamma*[1,1],[0,50],'color','red');hold off
    subplot(2,4,7);histogram(squeeze(posterSamples(nRound,:,3)));hold on;line(groundTruth.lapse*[1,1],[0,50],'color','red');hold off
    subplot(2,4,8);histogram(squeeze(posterSamples(nRound,:,4)));hold on;line(groundTruth.samples*[1,1],[0,50],'color','red');hold off
    pause()
    priorCZscore(i)= (groundTruth.priorC - mean(posterSamples(:,1))) / std(posterSamples(:,1));
    gammaZscore(i)= (groundTruth.gamma - mean(posterSamples(:,2))) / std(posterSamples(:,2));
    lapseZscore(i)= (groundTruth.lapse - mean(posterSamples(:,3))) / std(posterSamples(:,3));
    samplesZscore(i)= (groundTruth.samples - mean(posterSamples(:,4))) / std(posterSamples(:,4));
    
    % priorCError(i) = mean(posterSamples(:,1)) - groundTruth.priorC;%/ std(posterSamples(:,1));
    % gammaError(i) =  mean(posterSamples(:,2)) - groundTruth.gamma;% / std(posterSamples(:,2));
    % lapseError(i) =  mean(posterSamples(:,3)) - groundTruth.lapse;% / std(posterSamples(:,3));
    % samplesError(i) =  mean(posterSamples(:,4)) - groundTruth.samples;% / std(posterSamples(:,4));
    
end
figure
plot(priorCZscore,'-*');hold on;
plot(gammaZscore,'-*');hold on;
plot(lapseZscore,'-*');hold on;
plot(samplesZscore,'-*');hold on;
legend('priorCError','gammaError','lapseError','samplesError')
hold on
line([0,length(resNames)],[1.65,1.65],'linestyle','--','color','black')
line([0,length(resNames)],[-1.65,-1.65],'linestyle','--','color','black')
%%
%
[~,index] = sort(abs(priorCZscore)+abs(gammaZscore)+abs(lapseZscore),'ascend');

i = index(1);
name = [resFolder,'/',resNames(i).name];
fieldstofit = {'priorC','gamma','lapse','samples'};
% extract the ground truth of this dataset
for k = 1:numel(fieldstofit)
    ind =  strfind(name,fieldstofit{k}) + numel(fieldstofit{k});
    eval(sprintf('paraList = genDataInfo.%s_list;',fieldstofit{k}));
    eval(sprintf('groundTruth.%s = paraList(str2double(name(ind)));',fieldstofit{k}));
end
load(name);
figure
subplot(2,4,[1:4])
plot(log{1}.loss);
subplot(2,4,5);histogram(posterSamples(:,1));hold on;line(groundTruth.priorC*[1,1],[0,50],'color','red');hold off;xlim([0,1])
subplot(2,4,6);histogram(posterSamples(:,2));hold on;line(groundTruth.gamma*[1,1],[0,50],'color','red');hold off;xlim([0,1])
subplot(2,4,7);histogram(posterSamples(:,3));hold on;line(groundTruth.lapse*[1,1],[0,50],'color','red');hold off;xlim([0,1])
subplot(2,4,8);histogram(posterSamples(:,4));hold on;line(groundTruth.samples*[1,1],[0,50],'color','red');hold off;xlim([0,100])

i = index(end);
name = [resFolder,'/',resNames(i).name];
fieldstofit = {'priorC','gamma','lapse','samples'};
% extract the ground truth of this dataset
for k = 1:numel(fieldstofit)
    ind =  strfind(name,fieldstofit{k}) + numel(fieldstofit{k});
    eval(sprintf('paraList = genDataInfo.%s_list;',fieldstofit{k}));
    eval(sprintf('groundTruth.%s = paraList(str2double(name(ind)));',fieldstofit{k}));
end
load(name);
figure
subplot(2,4,[1:4])
plot(log{1}.loss);
subplot(2,4,5);histogram(posterSamples(:,1));hold on;line(groundTruth.priorC*[1,1],[0,50],'color','red');hold off;xlim([0,1])
subplot(2,4,6);histogram(posterSamples(:,2));hold on;line(groundTruth.gamma*[1,1],[0,50],'color','red');hold off;xlim([0,1])
subplot(2,4,7);histogram(posterSamples(:,3));hold on;line(groundTruth.lapse*[1,1],[0,50],'color','red');hold off;xlim([0,1])
subplot(2,4,8);histogram(posterSamples(:,4));hold on;line(groundTruth.samples*[1,1],[0,50],'color','red');hold off;xlim([0,100])