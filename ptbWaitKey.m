function [key, reaction_time, timeout] = ptbWaitKey(allowed_keys, timelimit)
if nargin < 2, timelimit = 0; end
timeout = false;
tstart = tic;
[~, ~, key] = KbCheck(-1);
while ~any(key == allowed_keys)
    [~, ~, key] = KbCheck(-1);
    if timelimit > 0 && toc(tstart) > timelimit
        timeout = true;
        break;
    end
end
reaction_time = toc(tstart);
key = intersect(find(key), allowed_keys);
end