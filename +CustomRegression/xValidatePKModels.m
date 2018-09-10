function [pkLL, regLL, expLL, linLL, reg_hpr] = xValidatePKModels(data, responses, folds)

[trials, frames] = size(data);

foldStart = round(linspace(1, trials+1, folds+1));
foldEnd = foldStart(2:end) - 1;

% Shuffle the data to avoid dependency of x-validation fold and time
shuffle = randperm(trials);
data = data(shuffle, :);
responses = responses(shuffle);

augData = [data ones(trials, 1)];

% Do regularization for ridge + AR2 separately so that same regularized hyperparams are used in
% every fold here.
reg_hpr = CustomRegression.xValidatePK(data, responses, logspace(-3, 3, 7), 0, logspace(-3, 3, 7), 1, 10);

for iFold=1:folds
    dataHoldOut = augData(foldStart(iFold):foldEnd(iFold), :);
    respHoldOut = responses(foldStart(iFold):foldEnd(iFold));
    dataFold = data([1:foldStart(iFold)-1 foldEnd(iFold)+1:end], :);
    respFold = responses([1:foldStart(iFold)-1 foldEnd(iFold)+1:end]);
    
    % Standard PK fit
    pk_w = CustomRegression.PsychophysicalKernel(dataFold, respFold, 0, 0, 0, 1);
    pkLL(iFold) = bernoulli_log_likelihood(dataHoldOut, respHoldOut, pk_w);
    
    % Regularized PK fit with a separate 10-fold XV for ridge + AR2 parameters
    reg_w = CustomRegression.PsychophysicalKernel(dataFold, respFold, reg_hpr(1), reg_hpr(2), reg_hpr(3), 1);
    regLL(iFold) = bernoulli_log_likelihood(dataHoldOut, respHoldOut, reg_w);
    
    % Exponential PK
    abb = CustomRegression.ExponentialPK(dataFold, respFold, 1);
    exp_weights = abb(1)*exp(abb(2)*(0:frames-1));
    linLL(iFold) = bernoulli_log_likelihood(dataHoldOut, respHoldOut, [exp_weights abb(3)]);
    
    % Linear PK
    sob = CustomRegression.LinearPK(dataFold, respFold, 1);
    lin_weights = sob(2) + (0:frames-1) * sob(1);
    expLL(iFold) = bernoulli_log_likelihood(dataHoldOut, respHoldOut, [lin_weights sob(3)]);
end

end


function LL = bernoulli_log_likelihood(data, responses, weights)
% Copied from CustomRegression.PsychophysicalKernel
logits = data * weights(:);
log_bernoulli = -responses(:) .* logits(:) + log(1 + exp(logits(:)));
LL = -sum(log_bernoulli);
end 