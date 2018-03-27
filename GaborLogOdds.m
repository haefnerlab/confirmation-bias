function [frame_odds, decision_odds] = GaborLogOdds(frames, left_category, right_category, contrast, noise, p_match)
%GABORLOGODDS comutes log odds that each frame is from the left or the
%right template, and the log odds that each frame is drawn from the left or
%right full-trial choice.
[n_imgs, h, w] = size(frames);
flat_frames = reshape(frames, [n_imgs, h*w]);

% TODO - rewrite for fourier domain signal.

diff_left = flat_frames - repmat(left_category(:)', n_imgs, 1);
diff_right = flat_frames - repmat(right_category(:)', n_imgs, 1);

log_left = -0.5 * sum(diff_left .* diff_left, 2) / noise;
log_right = -0.5 * sum(diff_right .* diff_right, 2)  / noise;

frame_odds = (log_left - log_right)';

% Derivation of log_decide_left:
%
% p_match * exp(log_left) + (1-p_match) * exp(log_right)
%   = exp(log_right) + p_match * [exp(log_left) - exp(log_right)]
% 
% note [exp(log_left) - exp(log_right)] = exp(log_sub(log_left, log_right))
% so, log(p_match * [exp(log_left) - exp(log_right)])
%  is log(p_match) + log_sub(log_left, log_right))
log_decide_left = real(log_add(log_right, log(p_match) + log_sub(log_left, log_right)));
log_decide_right = real(log_add(log_left, log(p_match) + log_sub(log_right, log_left)));

decision_odds = (log_decide_left - log_decide_right)';
end

function c = log_add(a, b)
% Computes c such that exp(c) = exp(a) + exp(b)
c = a + log(1 + exp(b - a));
c(a == -inf | exp(b - a) == inf) = b(a == -inf | exp(b - a) == inf);
c(b == -inf) = a(b == -inf);
end

function c = log_sub(a, b)
% Computes c such that exp(c) = exp(a) - exp(b)
c = a + log(1 - exp(b - a));
c(a == -inf | exp(b - a) == inf) = b(a == -inf | exp(b - a) == inf) + 1i*pi;
c(b == -inf) = a(b == -inf);
end