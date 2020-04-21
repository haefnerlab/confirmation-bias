function choice = psychoBernoulliSimulator(x,x_center,params)
% x is the input signal (stimulus) and should be a Nx1 vector
% x_center is the center value of stimulus
% params(1) should be sensitivity
% params(2) should be bias
% Another parameter is noise level of p. if this is not assigned by user,
% it would be set to default value
sensitivity = params(1);
bias = params(2);
if length(parms)<3
    % defalut value
    Noise_ratio = 0.05; 
end
if size(x,2) > 1
    error('x should be a N*1 vector');
end

N = size(x,1);
% probability of making choice as +1
p = (1 + exp(-  (bias+sensitivity*(x-x_center)) )).^-1;


p_noisy = p +  Noise_ratio * randn(N,1);
% boundary of p: should be within [0,1]
p_noisy(p_noisy<0) = 0;
p_noisy(p_noisy>1) = 1;

choice = zeros(N,1);
rnds = rand(N,1);
choice(rnds < p_noisy ) = 1;
