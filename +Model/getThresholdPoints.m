function [grid_pts, raw_pts] = getThresholdPoints(ps, params, pCorrect, nPoints)
%% Run model to get percent correct over full category/sensory information space

[ss, cc] = meshgrid(ps);
correct = Model.plotCategorySensorySpace(ps, ps, params, {}, [], false);

%% Drawing rays from the top right corner of the space, find equi-angle points along the threshold correct line

angles = linspace(0, pi/2, nPoints+2);
angles = angles(2:end-1);

t = linspace(0, 1);
for i=nPoints:-1:1
    s_ray = 1-t*cos(angles(i));
    c_ray = 1-t*sin(angles(i));
    
    correct_ray = interp2(ss, cc, correct, s_ray, c_ray);
    valid = ~isnan(correct_ray) & isfirstunique(correct_ray);

    best_t = interp1(correct_ray(valid), t(valid), pCorrect);
    
    raw_pts(i, :) = [1-best_t*cos(angles(i)) 1-best_t*sin(angles(i))];
end

%% Snap points to grid

for i=nPoints:-1:1
    [~, s_closest] = min(abs(ps - raw_pts(i, 1)));
    [~, c_closest] = min(abs(ps - raw_pts(i, 2)));
    
    grid_pts(i, :) = [ps(s_closest) ps(c_closest)];
end

end

function first = isfirstunique(values)
match = (values - values') == 0;
match(end+1,:) = true;
idxval = 1:length(values);
idxfirst = arrayfun(@(col) find(match(:,col), 1), idxval);
first = idxfirst == idxval;
first = reshape(first, size(values));
end