function rect = ptbCenteredRect(center, dims)
rect = [center(1) - dims(1) / 2, center(2) - dims(2) / 2, ...
    center(1) + dims(1) / 2, center(2) + dims(2) / 2];
rect = round(rect);
end

