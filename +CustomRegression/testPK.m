function testPK(trials, frames, ridge, ar1, curvature)

results_dir = fullfile('+CustomRegression', 'TestResults');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

regressors = randn(trials, frames);

    function compare(true_kernel, responses, name)
        [regression_weights, ~, errors, map_ridge, map_ar1, map_curvature] = CustomRegression.PsychophysicalKernel(regressors, responses, ridge, ar1, curvature);
        
        fig = figure(); hold on;
        plot(true_kernel, 'Color', 'k', 'LineWidth', 2);
        errorbar(regression_weights(1:end-1), errors(1:end-1));
        errorbar(frames+1, regression_weights(end), errors(end), 'Color', 'r', 'LineWidth', 2);
        title(name);
        saveas(fig, fullfile(results_dir, sprintf('%.2f %.2f %.2f %s.fig', map_ridge, map_ar1, map_curvature, name)));
    end

%% First test: flat kernel

kernel = ones(frames, 1);
responses = get_responses(kernel, regressors);
compare(kernel, responses, 'ideal');

%% Second test: decreasing PK

kernel = linspace(1, 0, frames)';
responses = get_responses(kernel, regressors);
compare(kernel, responses, 'linear decreasing');

%% Third test: uses only center values

kernel = zeros(frames, 1); kernel(end/2-1:end/2+1) = 1;
responses = get_responses(kernel, regressors);
compare(kernel, responses, 'middle three');

end

function responses = get_responses(kernel, regressors)
p = 1 ./ (1 + exp(-regressors * kernel));
responses = rand(size(p)) < p;
end