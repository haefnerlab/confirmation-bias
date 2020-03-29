function [M, L, U, median, weight_matrix] = BootstrapWeightsGabor(signals, choices, bootstrapsteps, hprs, normalize)
if nargin < 5, normalize = false; end

[trials, frames] = size(signals);
num_weights = frames + 1;
weight_matrix = zeros(bootstrapsteps, num_weights);

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi([1 trials], 1, trials);
    boot_choices = choices(index);
    boot_signals = signals(index, :);
   
    % Temporal PK regression
    weights = CustomRegression.PsychophysicalKernel(boot_signals, boot_choices, hprs(1), hprs(2), hprs(3), 1);
    if normalize
        weights(1:end-1) = weights(1:end-1) / mean(weights(1:end-1));
    end
    weight_matrix(i,:) = weights; 
end

[ M, L, U, median] = meanci(weight_matrix, .68);

end