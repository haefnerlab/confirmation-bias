function [key, reaction_time, timeout] = ptbWaitKey(allowed_keys, timelimit)
if nargin < 2, timelimit = 0; end
timeout = false;
tstart = tic;

[~, ~, key] = KbCheck(-1);
key = intersect(find(key), allowed_keys);
while isempty(key)
    [~, ~, key] = KbCheck(-1);
    key = intersect(find(key), allowed_keys);
    
    if timelimit > 0 && toc(tstart) > timelimit
        timeout = true;
        break;
    end
end
reaction_time = toc(tstart);
end