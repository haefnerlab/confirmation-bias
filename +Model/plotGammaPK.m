function [ fig ] = plotGammaPK(params, gammas)
%PLOTGAMMAPK Plot a PK for each of the given gamma values.

weights = cell(size(gammas));
errors = cell(size(gammas));

parfor i=1:length(gammas)
    params_copy = params;
    params_copy.gamma = gammas(i);
    
    [weights{i}, errors{i}, fig] = Model.plotPK(params_copy);
    close(fig);
end

%% Plot

fig = figure;
hold on;

colors = fadecolors([1 0 0], [0 0 1], length(gammas));

for i=1:length(gammas)
    errorbar(weights{i}(1:end-1), errors{i}(1:end-1), 'Color', colors(i,:), 'LineWidth', 2);
end
legend(arrayfun(@(g) ['\gamma = ' num2str(g)], gammas, 'UniformOutput', false), 'Location', 'best');
title('PK for different values of \gamma');

end