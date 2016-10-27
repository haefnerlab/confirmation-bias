function D = create_forward_difference_matrix(imsize, kernel)
%CREATE_FORWARD_DIFFERENCE_MATRIX constructs the forward finite-differnce
%matrix specified by the given kernel
%
% D = CREATE_FORWARD_DIFFERENCE_MATRIX(imsize, kernel) Constructs a matrix
% D such that, given an 'image' (or vector) I, 'reshape(D*I(:), size(I))'
% contains at each (x,y) location the given kernel applied to that
% location. For example, if I is [1 2 4; 7 11 12] and 'kernel' is [-1 1],
% then reshape(D*I(:), 2, 3) will be [1 2 0; 4 1 0].
%
% Kernel examples include [-1 1] for 1st (x) derivatives, [1 -2 1] for 2nd
% (x) derivatives or the transpose of these for y derivatives.

height = imsize(1);
width = imsize(2);
[kern_height, kern_width] = size(kernel);

tiles_x = width - kern_width + 1;
tiles_y = height - kern_height + 1;

% index + value pairs used to construct sparse matrix below
i_j_val = zeros(tiles_x*tiles_y, 3);

sub = 1;
for x=1:tiles_x
    for y=1:tiles_y
        % Computing the finite-difference estimate at (x, y)
        idx_xy = sub2ind(imsize, y, x);
        for kx=1:kern_width
            for ky=1:kern_height
                % Apply the kernel value at (kx, ky) to the 'pixel' at
                % (x+kx-1, y+ky-1)
                idx_kern = sub2ind(imsize, y+ky-1, x+kx-1);
                i_j_val(sub, :) = [idx_xy, idx_kern, kernel(ky, kx)];
                sub = sub + 1;
            end
        end
    end
end

D = sparse(i_j_val(:,1), i_j_val(:,2), i_j_val(:,3));

end