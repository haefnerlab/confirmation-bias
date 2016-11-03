function [Data, small_image_collection] = convertToSmallImageSize(Data, image_collection)

[trials, frames, ~, ~] = size(image_collection);

small_height = Data.image_length_y;
small_width  = Data.image_length_x;

%% Resize images in image_collection

small_image_collection = zeros(trials, frames, small_height, small_width);

for t=1:trials
    for f=1:frames
        small_image_collection(t, f, :, :) = ...
            imresize(squeeze(image_collection(t, f, :, :)), ...
            [small_height, small_width], 'nearest');
    end
end

%% Resize templates

if isfield(Data, 'left_template')
    Data.left_template  = imresize(Data.left_template,  [small_height, small_width], 'nearest');
    Data.right_template = imresize(Data.right_template, [small_height, small_width], 'nearest');
end

end