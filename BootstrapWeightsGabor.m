function [M, L, U] = BootstrapWeightsGabor(Test_Data, image_collection,automatic,bootstrapsteps)


weight_matrix=zeros(bootstrapsteps,Test_Data.image_length_y*Test_Data.image_length_x+Test_Data.number_of_images+1);

for i=1:bootstrapsteps
    disp(i)
    data={};
    data_image_collection=[];
    
    index=randi([1 Test_Data.current_trial],1,Test_Data.current_trial);
    data.move_on = Test_Data.move_on(index);
    data.step_size = Test_Data.step_size(index);
    data.reversal_counter = Test_Data.reversal_counter(index);
    data.contrast = Test_Data.contrast(index);
    data.ratio = Test_Data.ratio(index);
    data.correct_answer = Test_Data.correct_answer(index);
    data.staircase_answer = Test_Data.staircase_answer(index);
    data.reaction_time = Test_Data.reaction_time(index);
    data.choice = Test_Data.choice(index);
    data.accuracy = Test_Data.accuracy(index);
    data.order_of_orientations = Test_Data.order_of_orientations(index, :);
    data.log_odds = Test_Data.log_odds(index);
    data.average_orientations = Test_Data.average_orientations(:, index);
    data.image_template1 = Test_Data.image_template1(index, :);
    data.image_template2 = Test_Data.image_template2(index, :);
    data.image_template_difference = Test_Data.image_template_difference(index, :);
    if automatic==0
    data.eye_tracker_points = Test_Data.eye_tracker_points(index);
    end
    data.current_trial = Test_Data.current_trial;
    data.screen_frame = Test_Data.screen_frame;
    data.screen_resolution = Test_Data.screen_resolution;
    data.image_length_x = Test_Data.image_length_x;
    data.image_length_y = Test_Data.image_length_y;
    data.left_template = Test_Data.left_template;
    data.right_template = Test_Data.right_template;
    data.number_of_images = Test_Data.number_of_images;
    data_image_collection = image_collection(index, :, :, :);
    
    
    
    
    
    [trials, number_of_images, h, w] = size(data_image_collection);
    % number of trials, images shown per trial, and the height and width of the image in the experiment
    
    
    %========================================
    
    % Spatial + Temporal regression
    
    cell_images = mat2cell(data_image_collection, ones(trials, 1), number_of_images, h, w);
    cell_images = cellfun(@(im) permute(squeeze(im), [2 3 1]), cell_images, 'UniformOutput', false);
    [weights, ~, ~] = CustomRegression.PsychophysicalKernelImage(cell_images, data.choice, 0, 0, 10, 0, 0, 0);
    
    weight_matrix(i,:)=weights;
    
end

[ M, L, U ] = meanci( weight_matrix );

end