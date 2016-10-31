function image_array = makeImages(Data, order_of_orientations, contrast)

% makeImages is a function to determine what properties the gabor patches
% need to have before calling makeGabor to create the gabor patches. Then
% it passes the created gabors to trialStimuli to display the gabors on the
% screen

% order_of_orientations - an array of image orientations throughout the trial



res = Data.screen_resolution;

background = 127.0;

image_array = ones(Data.number_of_images, Data.image_length_x*res, Data.image_length_y*res);
for i = 1:Data.number_of_images
    if order_of_orientations(i) == 1
        image = Data.left_template * contrast + background;
    else
        image = Data.right_template * contrast + background;
    end
    static = [[ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)]];
   
    image = image + static;
    image(image > 255) = 255;
    image(image < 0) = 0;
    
    image_array(i,:,:) = image;
end

end