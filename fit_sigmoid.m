function [vm_params, ll] = fit_sigmoid(data,correct_incorrect)
%FIT_WEIBULL maximum likelihood fitting of a weibull cdf to the
%given data


% construct log likelihood function of params to be maximized
neg_log_likelihood = @(A) -log_likelihood(data,A(1),A(2),A(3));

% fminsearch to get params



ll = inf;
vm_params = [];

repeat=20;

for itr=1:repeat
   
    % reinitialize
    A0 = [rand,rand,rand];
    % solve
    [params, nll] = fminsearchbnd(neg_log_likelihood,A0,[0 -inf 0],[inf inf  1]); % max value for each VM parameter
    % check
    if nll < ll
        ll = nll;
        vm_params = params;
    end
    
end

function ll = log_likelihood(data,parameter1,parameter2,parameter3)
        
bernoulli_p = compute_sigmoid(data,parameter1,parameter2,parameter3);
ll = correct_incorrect.* log(bernoulli_p) + (1 - correct_incorrect) .* log(1 - bernoulli_p);
ll = sum(ll);
        
end

end