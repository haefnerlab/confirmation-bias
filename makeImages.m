function image_array = makeImages(Data)
%MAKEIMAGES creates noisy image frames for a single trial of the 'Gabor'
%experiment.
%
% images = MAKEIMAGES(Data) uses templates defined in the struct 'Data',
% using the left_tempate for each frame where
% order_of_orientations(current_trial,:) is 1 and the right_template
% otherwise. Returns 'images', a [frames x height x width] array of "low
% resolution" pixel values.
%
%  Data.image_length_{x,y}    - the width,height of the (low-res) images.
%  Data.{left,right}_template - the (size [y x]) tempates for the
%                               left/right categories.
%  Data.number_of_images      - frames per trial.

height = Data.image_length_y;
width = Data.image_length_x;
frames = Data.number_of_images;
t = Data.current_trial;
c = Data.contrast(t);
noise = Data.noise(t);

background = 127.0;
image_array = zeros(frames, height, width);
for i = 1:frames
    if Data.order_of_orientations(t, i) == 1
        image = Data.left_template * c + background;
    else
        image = Data.right_template * c + background;
    end
    
    % Add white pixel noise.
    image = image + noise * randn(height, width);
   
    % Clip pixel values to within the proper range.
    image(image > 255) = 255;
    image(image < 0) = 0;
    
    image_array(i,:,:) = image;
end
end