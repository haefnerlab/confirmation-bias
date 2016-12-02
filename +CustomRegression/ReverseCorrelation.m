function [ CTA, normed, smoothed ] = ReverseCorrelation(data, choices, smoothing)

if nargin < 3, smoothing = 5; end

stdev = sqrt(var(data));
CTA = mean(data(choices == 1, :)) - mean(data(choices == 0, :));
normed = CTA ./ stdev;
smoothed = smooth(normed, smoothing, 'moving');

end

