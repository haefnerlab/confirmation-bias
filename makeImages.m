function image_array = makeImages(Data, order_of_orientations, contrast)

% makeImages is a function to determine what porperties the gabor patches
% need to have before calling makeGabor to create the gabor patches. Then
% it passes the created gabors to trialStimuli to display the gabors on the
% screen

% order_of_orientations - an array of image orientations throughout the trial

%{
if Data.gabor_step > 1
    image_array = ones(Data.number_of_images, Data.gabor_step, Data.gabor_step);
else
    elements = floor((Data.gabor_endpoint - Data.gabor_startpoint) * (1/Data.gabor_step)+1);
    image_array = ones(Data.number_of_images, elements, elements);
end


for i = 1:Data.number_of_images
    gabor_angle = order_of_orientations(i);
    
    image = makeGabor(Data.gabor_sigma_x, Data.gabor_sigma_y, Data.gabor_spatial_frequency, Data.gabor_phase, ...
        gabor_angle, Data.gabor_startpoint, Data.gabor_step, Data.gabor_endpoint);
    
    if Data.gabor_step > 1
        image = (image .* contrast) + randn(Data.gabor_step, Data.gabor_step);
    else
        image = (image .* contrast) + randn(elements, elements);
    end
    
    image_array(i,:,:) = image;
end
%}

res = Data.screen_resolution;

background = 127.0;
contrast = 127.0 + contrast;

image_array = ones(Data.number_of_images, Data.image_length_x*res, Data.image_length_y*res);
for i = 1:Data.number_of_images
    if order_of_orientations(i) == 1
        image = [[ones(res,res)*background ones(res,res)*background ones(res,res)*contrast ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*contrast ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*contrast ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*contrast ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*contrast ones(res,res)*background ones(res,res)*background]];  % Left Orientation Binary Image
    else
        image = [[ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*contrast ones(res,res)*contrast ones(res,res)*contrast ones(res,res)*contrast ones(res,res)*contrast];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background];
                 [ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background ones(res,res)*background]];  % Right Orientation Binary Image
    end
    [x_axis, y_axis] = size(image);
    static = [[ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)];
             [ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0) ones(res,res)+(randn*16.0)]];
    %image = (image.*contrast) + 127;
    image = image + static;
    image(image > 255) = 255;
    image(image < 0) = 0;
    
    image_array(i,:,:) = image;
end

end