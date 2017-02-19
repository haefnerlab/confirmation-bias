function [cta, stderr] = CTA(data, responses)
%CTA compute difference of choice-triggered averages.

% Equalize variance for each regressor.
data = zscore(data);

% convert boolean to float type
u_responses = unique(responses);

assert(length(u_responses) == 2, 'Difference-of-CTAs only works with 2 choices');

data1 = data(responses == max(u_responses), :);
data2 = data(responses == min(u_responses), :);

cta = mean(data1, 1) - mean(data2, 1);
stderr = sqrt(var(data1, 1) + var(data2, 1)) / sqrt(size(data, 1) / 2);

end