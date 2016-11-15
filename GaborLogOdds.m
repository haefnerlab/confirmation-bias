function [frame_odds, decision_odds] = GaborLogOdds(frames, left_template, right_template, contrast, noise, p_match)
%GABORLOGODDS comutes log odds that each frame is from the left or the
%right template, and the log odds that each frame is drawn from the left or
%right full-trial choice.
t = GaborData.current_trial;

[n_imgs, h, w] = size(frames);
flat_frames = reshape(frames, [n_imgs, h*w]);

left_template = 127 + contrast * left_template;
right_template = 127 + contrast * right_template;

log_left = -0.5 * (flat_frames * left_template(:)) / noise;
log_right = -0.5 * (flat_frames * right_template(:)) / noise;

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
if a == -inf
    c = b;
elseif b == -inf
    c = a;
else
    c = a + log(1 + exp(b - a));
end
end

function c = log_sub(a, b)
% Computes c such that exp(c) = exp(a) - exp(b)
if a == -inf
    c = b;
elseif b == -inf
    c = a;
else
    c = a + log(1 - exp(b - a));
end
end