function showExampleImages(n_images, contrast, noise, p_match)
Data.image_length_y = 5;
Data.image_length_x = 5;
Data.number_of_images = n_images;
Data.left_template = eye(5);
Data.right_template = rot90(Data.left_template);
Data.pixel_noise = noise;

order = sign(p_match - rand(1:n_images));

images = makeImages(Data, order, contrast) / 256;

for i=1:n_images
    subplot(1, n_images, i);
    imagesc(squeeze(images(i, :, :)), [0 1]); colormap gray; axis image;
    xticks([]); yticks([]);
end
end