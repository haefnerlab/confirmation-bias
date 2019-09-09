function [M, L, U, median, weight_matrix, fits] = BootstrapExponentialWeightsGabor(Test_Data, bootstrapsteps, signalKappa, normalize)

if nargin < 3, signalKappa = 0; end
if nargin < 4, normalize = false; end

frame_signals = ComputeFrameSignals(Test_Data, signalKappa);

[trials, frames] = size(frame_signals);
weight_matrix = zeros(bootstrapsteps, frames);
biases = zeros(bootstrapsteps, 1);
fits = zeros(bootstrapsteps, 3);

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi([1 trials], 1, trials);
    boot_choices = Test_Data.choice(index) == +1;
    boot_signals = frame_signals(index, :);
   
    % Temporal PK regression
    abb = CustomRegression.ExponentialPK(boot_signals, boot_choices, 1);
    fits(i, :) = abb;
    weight_matrix(i, :) = abb(1)*exp(abb(2)*(0:frames-1));
    biases(i) = abb(3);
end

if normalize, weight_matrix = weight_matrix ./ mean(weight_matrix, 2); end
weight_matrix = [weight_matrix biases];
[ M, L, U, median] = meanci(weight_matrix, .68);

end