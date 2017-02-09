function fit_result = GaborPsychometric(Data, phase)
% Construct options for psignifit function.
options = struct;

if phase == 0
    options.sigmoidName  = 'weibull';
    options.expType      = '2AFC';
    
    % Count how often subject was correct at each contrast value.
    uniq_vals = unique(Data.contrast);
    yvals = arrayfun(@(u) mean(Data.accuracy(Data.contrast == u)), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(Data.contrast == u), uniq_vals);
elseif phase == 1
    options.sigmoidName  = 'norm';
    options.expType      = 'YesNo';
    
    % Count how often subject chose left at each ratio value.
    true_ratio = sum(Data.order_of_orientations, 2);
    uniq_vals = unique(true_ratio);
    
    % The next line exploits the fact that 'Left' is coded as 1 and 'right'
    % as 0, so mean() returns the fraction of left choices.
    yvals = arrayfun(@(u) mean(Data.choice(true_ratio == u)), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(true_ratio == u), uniq_vals);
end

% Run PM fitting.
fit_result = psignifit([uniq_vals(:) yvals(:) num_trials_at_vals(:)], options);
end