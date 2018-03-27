function [M, L, U, median, weight_matrix] = BootstrapCTA(Test_Data, bootstrapsteps, signalKappa)

if nargin < 3, signalKappa = 0; end

frame_signals = ComputeFrameSignals(Test_Data, signalKappa);

[trials, frames] = size(frame_signals);
weight_matrix = zeros(bootstrapsteps, frames);

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi([1 trials], 1, trials);
    boot_choices = Test_Data.choice(index) == +1;
    boot_signals = frame_signals(index, :);
   
    % Compute CTA
    weight_matrix(i,:) = mean(boot_signals(boot_choices, :)) - ...
        mean(boot_signals(~boot_choices, :)); 
end

% Force 'bias' to 0
% TODO
weight_matrix(:, end+1) = 0;

[ M, L, U, median ] = meanci(weight_matrix, .68);

end