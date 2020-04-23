function pairs = GridSpace(alpha,beta,sigma)
% function to generate grid points in parameter space
% (should be saved in the samefolder with analytical_posterior.m)

% Input: 
    % alpha: m1-by-1, 
    % beta: m2-by-1, 
    % sigma: m3-by-1
% Output: 
    % params_space: M-by-1 cells, 
    % each element: {1,4,20}--> {1*3} represent a specific pair of {alpha,beta,sigma}.
    
    % pairs: M-by-1 matrix, 
    
[param_x,param_y,param_z] = meshgrid(alpha,beta,sigma);
pairs = [param_x(:),param_y(:),param_z(:)]; % demension: M-by-1, where M= (m1)*(m2)*(m3)
% params_space = num2cell(pairs,2); % changing pairs mat [1,3] --> cell {1,3}
end

