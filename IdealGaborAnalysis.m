function [] = IdealGaborAnalysis(subjectID,groupings,directory)




prelimFile = [directory 'RawData/' subjectID '-GaborDataRatio.mat'];
if ~exist(prelimFile, 'file')
    disp(strcat('ERROR! Missing File: ', prelimFile));  % Return an error message for missing file
    disp(strcat('Maybe the Preliminary phase is saved under a different name?'));
    return;
else
    results=load(prelimFile); % Load Preliminary_Data
end
test_collection_of_images = results.image_collection;

test_choice = results.Preliminary_Data.choice;
[trials, number_of_images, h, w] = size(test_collection_of_images);
% number of trials, images shown per trial, and the height and width of the image in the experiment

sublength = number_of_images / groupings;
%% Without bootstrap
cell_images = mat2cell(test_collection_of_images, ones(trials, 1), number_of_images, h, w);
cell_images = cellfun(@(im) permute(squeeze(im), [2 3 1]), cell_images, 'UniformOutput', false);
[weights, ~, errors] = ...
    CustomRegression.PsychophysicalKernelImage(cell_images, test_choice, 0, 0, 10, 0, 0, 0);

Get_Figure('Temporal Kernel');hold on;
errorbar(1:number_of_images, weights(h*w+1:end-1), errors(h*w+1:end-1));  % Blue plot
plot([mean(reshape(1:number_of_images, [sublength groupings]))], [sum(reshape(weights(h*w+1:end-1), [sublength groupings]))]/sublength,'r*-');    % Red plot


title('Weighting the Image Frames')
legend('Actual Weight in each frame','Actual weights summed in groups','Location','southoutside')




Get_Figure('Image Kernel');
r=reshape(weights(1:h*w), [h w]);
imagesc(r);
colorbar;
end