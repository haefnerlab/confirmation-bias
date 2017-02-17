function [ fig ] = plotGammaPK(trials, frames, params, gammas)
%PLOTGAMMAPK Plot a PK for each of the given gamma values.

% All model versions can use the exact same data.
[data, data_prefix] = Model.genDataWithParams(trials, frames, params);

weights = cell(size(gammas));
errors = cell(size(gammas));

parfor i=1:length(gammas)
    params_copy = params;
    params_copy.gamma = gammas(i);
    
    results = Model.loadOrRunSamplingModel(data, data_prefix, params_copy);
    
    [w, e] = Model.loadOrRunModelPK(...
        Model.getModelStringID(data_prefix, params_copy), ...
        data, results, [1 0 10]);
    
    weights{i} = w;
    errors{i} = e;
end

%% Plot

fig = figure;
hold on;

colors = winter(length(gammas));

for i=1:length(gammas)
    errorbar(weights{i}, errors{i}, 'Color', colors(i,:), 'LineWidth', 2);
end
legend(arrayfun(@(g) ['\gamma = ' num2str(g)], gammas, 'UniformOutput', false));
title('PK for different values of \gamma');

end

