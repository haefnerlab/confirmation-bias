function image_array = makeImages(Data, order_of_orientations, contrast)
%MAKEIMAGES creates noisy image frames for a single trial of the 'Gabor'
%experiment.
%
% images = MAKEIMAGES(Data, order_of_orientations, contrast) uses templates
% defined in the struct 'Data', using the left_tempate for each frame where
% order_of_orientations is 1 and the right_template otherwise. Returns
% 'images', a [frames x height x width] array of "low resolution" pixel
% values.
%
%  Data.image_length_{x,y}    - the width,height of the (low-res) images.
%  Data.{left,right}_template - the (size [y x]) tempates for the
%                               left/right categories.
%  Data.number_of_images      - frames per trial.

height = Data.image_length_y;
width = Data.image_length_x;
frames = Data.number_of_images;

background = 127.0;
image_array = zeros(frames, height, width);
for i = 1:frames
    if order_of_orientations(i) == 1
        image = Data.left_template * contrast + background;
    else
        image = Data.right_template * contrast + background;
    end
    
    % Add white pixel noise.
    image = image + Data.pixel_noise * randn(height, width);
   
    % Clip pixel values to within the proper range.
    image(image > 255) = 255;
    image(image < 0) = 0;
    
    image_array(i,:,:) = image;
end
end