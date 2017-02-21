function [val, idx] = closest(arr, val)
%CLOSEST Find the closest value and its index in an array.
diffs = abs(arr - val);
[~, idx] = min(diffs);
val = arr(idx);
end