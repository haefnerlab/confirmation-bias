function [M, L, U, median, fits_matrix] = BootstrapExponentialWeightsGabor(Test_Data, bootstrapsteps, signalKappa, binarize)

if nargin < 3, signalKappa = 0; end
if nargin < 4, binarize = false; end

frame_signals = ComputeFrameSignals(Test_Data, signalKappa);
if binarize
    frame_signals = sign(frame_signals);
end

[trials, frames] = size(frame_signals);
num_weights = frames + 1;
fits_matrix = zeros(bootstrapsteps, 2);

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi([1 trials], 1, trials);
    boot_choices = Test_Data.choice(index) == +1;
    boot_signals = frame_signals(index, :);
   
    % Temporal PK regression
    weights = CustomRegression.PsychophysicalKernel(boot_signals, boot_choices, 0, 0, 0, 1);
    fits_matrix(i,:) = CustomRegression.expFit(weights);
end

[ M, L, U, median] = meanci(fits_matrix, .68);

end