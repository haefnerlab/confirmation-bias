function key = ptbWaitKey(allowed_keys)
[~, ~, key] = KbCheck(-1);
while ~any(key == allowed_keys)
    [~, ~, key] = KbCheck(-1);
end
end