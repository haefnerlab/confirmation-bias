function [M, L, U, weight_matrix] = BootstrapWeightsGabor(Test_Data, bootstrapsteps, signalKappa, binarize)

if nargin < 3, signalKappa = 0; end
if nargin < 4, binarize = false; end 

memodir = fullfile(pwd, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

frame_signals = ComputeFrameSignals(Test_Data, signalKappa);
if binarize
    frame_signals = sign(frame_signals);
end

[trials, frames] = size(frame_signals);
num_weights = frames + 1;
weight_matrix = zeros(bootstrapsteps, num_weights);

parfor i=1:bootstrapsteps
    % Randomly resample trials with replacement
    index = randi([1 trials], 1, trials);
    boot_choices = Test_Data.choice(index) == +1;
    boot_signals = frame_signals(index, :);
   
    % Temporal PK regression
    weights = CustomRegression.PsychophysicalKernel(boot_signals, boot_choices, 1, 0, 10, false, zeros(1, frames));
    weight_matrix(i,:) = weights; 
end

[ M, L, U ] = meanci(weight_matrix, .68);

end