function [M, L, U, weight_matrix] = BootstrapLinearPKFit(Test_Data, bootstrapsteps, signalKappa, binarize)

if nargin < 3, signalKappa = 0; end
if nargin < 4, binarize = false; end 

frame_signals = ComputeFrameSignals(Test_Data, signalKappa);
if binarize
    frame_signals = sign(frame_signals);
end

[trials, frames] = size(frame_signals);
weight_matrix = zeros(bootstrapsteps, frames);
biases = zeros(bootstrapsteps, 1);

for i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi(trials, 1, trials);
    boot_choices = Test_Data.choice(index) == +1;
    boot_signals = frame_signals(index, :);
   
    % Temporal PK regression
    sob = CustomRegression.LinearPK(boot_signals, boot_choices);
    weight_matrix(i, :) = sob(2) + (0:frames-1) * sob(1); 
    biases(i) = sob(3);
end
weight_matrix = [weight_matrix biases];
[ M, L, U ] = meanci(weight_matrix, .68);

end