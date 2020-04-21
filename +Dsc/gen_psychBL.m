function y = gen_psychBL(x,x_mean, params)
% fuction to generate bonulli response (0/1 binary choice)
% INPUT
% x: 1*M vector input signals
% N: number of samples to generate
% parms: [sensitivity,bias,Noise_ratio(default = 0.05)]
% first two terms are important to define psychmetric curve
?
?
%% Step 1 data generation
% 1.1 sigmod function
% 1.2 add some random noise
% 1.3 bin the response and get the curve
?
% ---Make up some data
% x= 1:0.1:10;
% sensitivity index: higher more sensitive
% sensitivity = 1;
% bias
% bias= 5;
% ---
?
sensitivity = params(1);
bias = params(2);
?
if length(parms)<3
    % defalut value
    Noise_ratio = 0.05; 
end
?
%psych_fun = @(x) (1 + exp(- sensitivity * (x - bias))).^-1;
?psych_fun = @(x) (1 + exp(- (bias+sensitivity * (x - x_mean)))).^-1;
% rand sample from the curve
rng('default')
sample_index = randsample(length(x),N,'true');
sample_x = x(sample_index);
?
% control the Noise_ratio: how noisy the samples are (larger more noisy)
rng('default')
% Noise_ratio = 0.05; 
sample_y = psych_fun(sample_x) + Noise_ratio * normrnd(0,1,[size(sample_x,1),size(sample_x,2)]);
?
% set the boundary of the sample: [0,1]
sample_y( find(sample_y>1) ) = 1; sample_y(find(sample_y<0)) = 0; 
?
%% Step 2 generate the response based on Bernoulli(p): p = sample_y
y = binornd(1,sample_y);
