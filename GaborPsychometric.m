function [fit_result, uniq_vals, yvals, stderrs] = GaborPsychometric(Data, phase)
% Construct options for psignifit function.
options = struct;

if phase == 0
    options.sigmoidName  = 'weibull';
    options.expType      = '2AFC';
    
    % Count how often subject was correct at each contrast value.
    uniq_vals = unique(Data.contrast);
    yvals = arrayfun(@(u) mean(Data.accuracy(Data.contrast == u)), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(Data.contrast == u), uniq_vals);
    stderrs = arrayfun(@(u) std(Data.accuracy(Data.contrast == u)), uniq_vals) ./ sqrt(num_trials_at_vals);
elseif phase == 1
    options.sigmoidName  = 'norm';
    options.expType      = 'YesNo';
    
    % Count how often subject chose left at each ratio value.
    uniq_vals = unique(Data.true_ratio);
    
    % The next line exploits the fact that 'Left' is coded as 1 and 'right'
    % as 0, so mean() returns the fraction of left choices.
    yvals = arrayfun(@(u) mean(Data.choice(Data.true_ratio == u)), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(Data.true_ratio == u), uniq_vals);
    stderrs = arrayfun(@(u) std(Data.choice(Data.true_ratio == u)), uniq_vals) ./ sqrt(num_trials_at_vals);
elseif phase == 2
    options.sigmoidName  = 'norm';
    options.expType      = '2AFC';
    
    % Count how often subject was correct at each kappa value.
    uniq_vals = unique(Data.noise);
    yvals = arrayfun(@(u) mean(Data.accuracy(Data.noise == u)), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(Data.noise == u), uniq_vals);
    stderrs = arrayfun(@(u) std(Data.accuracy(Data.noise == u)), uniq_vals) ./ sqrt(num_trials_at_vals); 
end

% Run PM fitting.
fit_result = psignifit([uniq_vals(:) yvals(:) num_trials_at_vals(:)], options);
end