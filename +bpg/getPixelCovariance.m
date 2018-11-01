function [C, K] = getPixelCovariance(width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix)
%BPG.GETPIXELCOVARIANCE given parameters of BPG stimulus, compute the pixel-wise covariance matrix
%such that BPG images could be drawn from mvnrnd(zeros(...), C)
%
%C = BPG.GETPIXELCOVARIANCE(width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix) computes
%[width*width x width*width] covariance matrix C.
%
% [C, K] = BPG.GETPIXELCOVARIANCE(...) also returns K such that C=K'*K. Using K can be more
% efficient than using C in some cases.

px = width * width;

[I, ~, filterF, aperture] = bpg.genImages(1, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix);

filterPx = fftshift(real(ifft2(ifftshift(filterF))));
filterPx = filterPx(:);

K = eye(px); 
for i=1:px
    K(:,i) = filterPx([i:end 1:i-1]);
end

K = K .* aperture(:);
C = K * K';

% figure;
% r = [-1.8e-3 1.8e-3];
% subplot(1, 3, 1); imagesc(squeeze(I), r); axis image; colormap gray; colorbar;
% subplot(1, 3, 2); imagesc(reshape(K * randn(px, 1), [width width]), r); axis image; colormap gray; colorbar;
% subplot(1, 3, 3); imagesc(reshape(mvnrnd(zeros(px, 1), C), [width width]), r); axis image; colormap gray; colorbar;

end