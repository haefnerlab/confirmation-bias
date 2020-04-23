
%% test
% clear
% x = 20*[1:10]
% y = 2*x+3+randn(size(x))
% data = [x',y'];
% param_domin = [0,10; 0,10; 0,1];
% binsz = 0.1;
% [a,b,s]= analytical_posterior(data,param_domin,binsz);
% subplot(1,3,1)
% plot(a)
% subplot(1,3,2)
% plot(b)
% subplot(1,3,3)
% plot(s)
%%
function [post_a,post_b,post_s] = analytical_posterior(data,param_domin,binsz)
% Input:
% data: [x,y] should be n-by-1 vector
% x = [x1,...,xn]', y = [y1,...,yn]'
% prior: m-by-3 vector

% Output: (m-by-3 vector)
% posterior: for all the point in the parameter space
xx = data(:,1);
yy = data(:,2);

% param: [alpha,beta,sigma] should be m-by-1 vector
alpha_space = linspace (param_domin(1,1),param_domin(1,2), abs(param_domin(1,1)- param_domin(1,2))/binsz);
beta_space = linspace (param_domin(2,1),param_domin(2,2), abs(param_domin(2,1)- param_domin(2,2))/binsz);
SIG_space = linspace (param_domin(3,1),param_domin(3,2), abs(param_domin(3,1)- param_domin(3,2))/binsz);
% sanity check: sigma can't be negative (lower bound is 0)
SIG_space ( (SIG_space<0) ) = 0;

% grids in parms space
pairs = GridSpace(alpha_space,beta_space,SIG_space);

% likelihood function:
% vary the pair of params {alpha_0,beta_0,sigma_0} in params space
for i = 1: length(pairs)
    alpha = pairs(i,1);
    beta = pairs(i,2);
    SIG = pairs(i,3);
    % for a specific (fixed) pair of params: {alpha_0,beta_0,sigma_0}
    % compute the likelihood: for [x1,...,xn]', [y1,...,yn]'
     mu = (alpha * xx + beta);
    n_joint_likeli(:,i) = normpdf(yy,mu,SIG); 
end

% replace the NaN value with 0
n_joint_likeli(isnan(n_joint_likeli)) = 0;


% need cumprod: \pi p(xn,yn| alpha_i,beta_i,sigma_i): 
% for a sepcific params {alpha_i,beta_i,sigma_i}
before_cumprod_joint = cumprod(n_joint_likeli);
joint_likelihood = before_cumprod_joint (end,:)';

% likelihood for alpha beta sigma
[uval,~,subs] = unique(pairs(:,1));
LL_a = accumarray(subs,joint_likelihood);

[uval,~,subs] = unique(pairs(:,2));
LL_b = accumarray(subs,joint_likelihood);

[uval,~,subs] = unique(pairs(:,3));
LL_s = accumarray(subs,joint_likelihood);

% normalized posterior
post_a = LL_a / sum(LL_a);
post_b = LL_b / sum(LL_b);
post_s = LL_s / sum(LL_s);
end









