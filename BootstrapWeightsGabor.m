function [M, L, U] = BootstrapWeightsGabor(Test_Data, image_collection, bootstrapsteps, ideal_template)

[trials, frames, h, w] = size(image_collection);
assert(Test_Data.current_trial == trials, 'Mismatch between Test_Data.current_trial and size(image_collection, 1)');
num_weights = h*w + frames + 1;
weight_matrix = zeros(bootstrapsteps, num_weights);

template_difference = Test_Data.left_template(:) - Test_Data.right_template(:);

% Permute dimensions of images to shape [trials h w frames] once so it
% doesn't have to be done each loop.
if ~ideal_template
    image_collection = permute(image_collection, [1 3 4 2]);
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
        img_rows = reshape(boot_images, [trials*frames, h*w]);
        img_data = reshape(img_rows * template_difference, [trials, frames]);
        weights = vertcat(template_difference, ...
            CustomRegression.PsychophysicalKernel(img_data, boot_choices, 1, 0, 10));
    else
        cell_images = mat2cell(boot_images, ones(trials, 1), h, w, frames);
        cell_images = cellfun(@squeeze, cell_images, 'UniformOutput', false);
        weights = CustomRegression.PsychophysicalKernelImage(cell_images, boot_choices, ...
            10, 0, 10, 0, 0, 0, template_difference);
    end
    weight_matrix(i,:) = weights; 
end

[ M, L, U ] = meanci( weight_matrix );

end