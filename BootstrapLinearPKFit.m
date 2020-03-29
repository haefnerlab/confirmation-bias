function [M, L, U, median, weight_matrix, fits] = BootstrapLinearPKFit(signals, choices, bootstrapsteps, normalize)

if nargin < 4, normalize = false; end 

[trials, frames] = size(signals);
weight_matrix = zeros(bootstrapsteps, frames);
biases = zeros(bootstrapsteps, 1);
fits = zeros(bootstrapsteps, 3);

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi(trials, 1, trials);
    boot_choices = choices(index) == +1;
    boot_signals = signals(index, :);
   
    % Temporal PK regression
    sob = CustomRegression.LinearPK(boot_signals, boot_choices, 1);
    fits(i, :) = sob;
    weight_matrix(i, :) = sob(2) + (0:frames-1) * sob(1); 
    biases(i) = sob(3);
end
if normalize
    % Normalizing means scaling down both the offset and the slope terms
    fits(:,1) = fits(:,1) ./ mean(weight_matrix, 2);
    fits(:,2) = fits(:,2) ./ mean(weight_matrix, 2);
    weight_matrix = weight_matrix ./ mean(weight_matrix, 2);
end
weight_matrix = [weight_matrix biases];
[ M, L, U, median] = meanci(weight_matrix, .68);

end