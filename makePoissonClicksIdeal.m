function [number_of_left_clicks, number_of_right_clicks, y,l,r] = makePoissonClicksIdeal(bins, left_click_rate, right_click_rate, leftVolume, rightVolume)

%sampling_rate = 5000;    % Sampling rate for playing audio
%left_click_rate = 10;    % Clicks per sec
%right_click_rate = 60;   % Clicks per sec
%stimulus_duration = 1; % sec

lambda_left = left_click_rate/bins; % # clicks per bin
lambda_right = right_click_rate/bins; % # clicks per bin


clickTrainLeft = double(binornd(1, lambda_left, 1, bins));
clickTrainRight = double(binornd(1, lambda_right, 1, bins));
l=clickTrainLeft;
r=clickTrainRight;
number_of_left_clicks = sum(clickTrainLeft);     % Tells us how many clicks for the left and right ear
number_of_right_clicks = sum(clickTrainRight);

clickTrainLeft = clickTrainLeft .* leftVolume;    % Determines how loud to make the clicks
clickTrainRight = clickTrainRight .* rightVolume;

y = [clickTrainLeft; clickTrainRight];
end


% Each 'click' was a sum of pure tones (at 2, 4, 6, 8, and 16 kHz) convolved with a cosine envelope 3 msec in width
% Click rate is L + R = 40 clicks/sec
% Options - 39:1, 37:3, 31:9, or 26:14 ratio
% Memory Gap? Pause between end of trial and before answering

% Log base e of 39/1 = 3.6635616461
% Log base e of 37/3 = 2.512305624
% Log base e of 31/9 = 1.2367626271
% Log base e of 26/14 = 0.61903920841