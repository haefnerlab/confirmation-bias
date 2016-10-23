function testPK(trials, ridge, ar1)

left = 0.1 * round(rand(trials, 120));
right = 0.1 * round(rand(trials, 120));

regressors = [left right ones(trials, 1)];

    function compare(true_left, true_right, responses, name)
        [regression_weights, ~, errors] = CustomRegression.PsychophysicalKernel(regressors, responses, ridge, ar1, true);
        
        figure(); hold on;
        est_left = regression_weights(1:120);
        est_right = regression_weights(121:end-1);
        plot(true_left);
        plot(-true_right);
        errorbar(est_left, errors(1:120));
        errorbar(est_right, errors(121:240));
        title(name);
    end

%% First test: flat both

kernel = ones(120, 1);
responses = left * kernel - right * kernel > 0;
compare(kernel, kernel, responses, 'ideal');

%% Second test: flat one side, other unused, but a bias

kernel = ones(120, 1);
responses = left * kernel > 6; % 120 / 2 * 0.1
compare(kernel, zeros(size(kernel)), responses, 'left only');

%% Third test: decreasing PK

kernel = linspace(1, 0, 120)';
responses = left * kernel - right * kernel > 0;
compare(kernel, kernel, responses, 'linear decreasing');

%% Fourth test: uses only center values

kernel = zeros(120, 1); kernel(end/2-1:end/2+1) = 1;
responses = left * kernel - right * kernel > 0;
compare(kernel, kernel, responses, 'middle three');

end