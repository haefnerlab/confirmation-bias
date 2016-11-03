function small_contents = convertToSmallImageSize(old_file, new_file, data_field, image_field)

contents = load(old_file);
small_contents = contents;

[trials, frames, ~, ~] = size(contents.(image_field));

small_height = contents.(data_field).image_length_y;
small_width  = contents.(data_field).image_length_x;

%% Resize images in contents.(image_field)
small_contents.(image_field) = zeros(trials, frames, small_height, small_width);

for t=1:trials
    for f=1:frames
        small_contents.(image_field)(t, f, :, :) = ...
            imresize(squeeze(contents.(image_field)(t, f, :, :)), ...
            [small_height, small_width], 'nearest');
    end
end

%% Resize templates
if isfield(contents.(data_field), 'left_template')
    small_contents.(data_field).left_template  = imresize(contents.(data_field).left_template,  [small_height, small_width], 'nearest');
    small_contents.(data_field).right_template = imresize(contents.(data_field).right_template, [small_height, small_width], 'nearest');
end

%% Save converted data to the new file
save(new_file, '-struct', 'small_contents');

end