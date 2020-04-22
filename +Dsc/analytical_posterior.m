rng('default')
kk= 20;
param = rand(kk,3);

xx = 10*[1:100]';
yy = (xx * param(1,1) + param(1,2)) +rand(1);

data = [xx, yy];

prior = ones(kk,1)/sum(ones(kk,1));

% function posterior = analytical_posterior(data,param,prior)
% Input:
    % data: [x,y] should be n-by-1 vector
    % x = [x1,...,xn]', y = [y1,...,yn]'
    % param: [alpha,beta,sigma] should be m-by-1 vector
    % prior: m-by-1 vector
% Output: (m-by-1 vector)
    % posterior: for all the point in the parameter space
x = data(:,1);
y = data(:,2);
alpha = param(:,1);
beta = param(:,2);
sigma = param(:,3);
% length of the params
Len = length(alpha);
% initial value before doing product
norm_likelihood = ones(1,Len); % 1-by-m vector
for i = 1: Len  
    for n = 1: length(y)
        mu = alpha(i) * x(n) + beta(i);
        % tmp = p(xn,yn| theta_i)
        tmp = exp(-0.5 * ((y(n) - mu)/sigma(i)).^2) / (sqrt(2*pi) * sigma(i));
        % likelihood = p(x1,y1| theta_i) * p(x2,y2| theta_i) *... p(xn,yn| theta_i)
        norm_likelihood(i) = norm_likelihood(i) * tmp;
    end
end
% posterior: m-by-1 vector
unnormalized_post = norm_likelihood'.* prior;
posterior = unnormalized_post /sum(unnormalized_post);
% test
