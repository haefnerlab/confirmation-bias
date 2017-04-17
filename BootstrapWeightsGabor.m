function [M, L, U] = BootstrapWeightsGabor(Test_Data, bootstrapsteps, ideal_template)

[trials, frames] = size(Test_Data.frame_categories);
sz = Test_Data.stim_size;
num_weights = sz*sz + frames + 1;
weight_matrix = zeros(bootstrapsteps, num_weights);

% TODO - rewrite for 1-param fourier domain spatial kernel.
template_difference = Test_Data.left_template(:) - Test_Data.right_template(:);

% Regenerate stimulus frames.
% TODO - smarter use of space for high res images.
image_collection = zeros(trials, sz, sz, frames);
parfor t=1:trials
    image_collection(t, :) = permute(GaborStimulus(Test_Data, t), [2 3 1]);
end

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index=randi([1 trials], 1, trials);
    boot_choices = Test_Data.choice(index) == +1;
    boot_images = image_collection(index, :, :, :);
   
    % Spatial + Temporal regression
    if ideal_template
        % Convolve with the ideal template to get 1d signal over time for
        % each trial.
        img_rows = reshape(boot_images, [trials*frames, sz*sz]);
        img_data = reshape(img_rows * template_difference, [trials, frames]);
        weights = vertcat(template_difference, ...
            CustomRegression.PsychophysicalKernel(img_data, boot_choices, 1, 0, 10));
    else
        cell_images = mat2cell(boot_images, ones(trials, 1), sz, sz, frames);
        cell_images = cellfun(@squeeze, cell_images, 'UniformOutput', false);
        weights = CustomRegression.PsychophysicalKernelImage(cell_images, boot_choices, ...
            1, 0, 0, 0, 0, 0, template_difference);
    end
    weight_matrix(i,:) = weights; 
end

[ M, L, U ] = meanci(weight_matrix, .68);

end