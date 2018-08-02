%% Demo showing how BPG stimulus is MVN-distributed in pixel space

frames = 4;
sz = 100;
mu_rho = 0.1194;
std_rho = 0.0597;
kappa = 0.8;
annulus = 25;

% Note: to show that scales are comparable, need to comment out 'im = im / max(abs(im(:)));' line in
% bpg.genImages
[im, ~, filterF, aperture] = bpg.genImages(frames, sz, mu_rho, std_rho, 30, kappa, annulus);
range = max(abs(im(:)));

subplot(2,frames+1,1);
imagesc(abs(filterF));
axis image; set(gca, 'YDir', 'normal'); axis off;
title('Fourier domain mask');

for f=1:frames
    subplot(2,frames+1,f+1);
    imagesc(squeeze(im(f,:,:)), [-range range]);
    axis image; set(gca, 'YDir', 'normal'); axis off;
    colorbar;
end

%% Above: fourier version. Below: MVN version

filterPx = real(fftshift(ifft2(ifftshift(filterF))));

subplot(2,frames+1,frames+2);
imagesc(filterPx);
axis image; set(gca, 'YDir', 'normal'); axis off;
title('F[mask] in pixel-space');

px = sz^2;
kernel = zeros(px, px);
for i=1:sz
    for j=1:sz
        % Get 'impulse response' of filterPx at location (i,j)
        inpt = zeros(sz, sz);
        inpt(i,j) = 1;
        resp = conv2(inpt, filterPx, 'same');
        kernel(:, sub2ind([sz sz], i, j)) = resp(:);
    end
end

% NOTE: pixel covariance is kernel' * kernel. We could compute this and pass it to mvnrnd, but then
% mvnrnd would just decompose Cov back into the kernel. We can instead pass the cholesky
% decomposition of the covariance matrix directly into mvnrnd as an undocumented 4th argument while
% passing a 'dummy' matrix in for the real covariance.
covChol = aperture(:) .* kernel;
dummy_cov = speye(px, px);
im = mvnrnd(zeros(px, 1), dummy_cov, frames, covChol);
range = max(abs(im(:)));
for f=1:frames
    subplot(2,frames+1,frames+2+f);
    imagesc(reshape(im(f, :), [sz sz]), [-range range]);
    axis image; set(gca, 'YDir', 'normal'); axis off;
    colorbar;
end