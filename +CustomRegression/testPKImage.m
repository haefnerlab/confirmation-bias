function testPKImage(trials, w, frames, ridge, ar1, curvature, ridge_sp, ar1_sp, curvature_sp)

results_dir = fullfile('+CustomRegression', 'TestResults');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

trialdata = arrayfun(@(t) randn(w, w, frames), 1:trials, 'UniformOutput', false);
flatdata = cellfun(@(trial) reshape(trial, w*w, frames), trialdata, 'UniformOutput', false);

    function compare(template_left, template_right, temporal, bias, name)
        % simulate results
        diff = template_left(:)-template_right(:);
        timedata = cellfun(@(trial) diff'*trial, flatdata, 'UniformOutput', false);
        timedata = vertcat(timedata{:});
        responses = timedata*temporal + bias > 0;

        % run regression
        [regression_weights, ~, errors] = CustomRegression.PsychophysicalKernelImage(trialdata, responses, ridge, ar1, curvature, ridge_sp, ar1_sp, curvature_sp);
        
        fig = figure();
        % plot spatial weights
        subplot(1,3,1); imagesc(template_left-template_right); colormap gray; axis image; colorbar;
        subplot(1,3,2); imagesc(reshape(regression_weights(1:w*w),w,w)); colormap gray; axis image; colorbar;
        % plot temporal weights
        subplot(1,3,3); hold on;
        plot(temporal);
        if mean(regression_weights(w*w+1:end-1)) > 0
            errorbar(regression_weights(w*w+1:end-1), errors(w*w+1:end-1));
        else
            errorbar(-regression_weights(w*w+1:end-1), errors(w*w+1:end-1));
        end
            ylim([0 inf]);
        % label and save
        title(name);
        saveas(fig, fullfile(results_dir, sprintf('Image - %.2f %.2f %.2f %.2f %.2f %.2f %s.fig', ridge, ar1, curvature, ridge_sp, ar1_sp, curvature_sp, name)));
    end

%% First test: bar templates with flat temporal

temporal = ones(frames, 1);
midpt = round(w/2);
template_left = zeros(w); template_left(midpt-1:midpt+1, :) = 1;
template_right = template_left';

compare(template_left, template_right, temporal, 0, 'bars & flat');

%% Second test: bar templates with decreasing temporal

temporal = linspace(1,0,frames)';
template_left = zeros(w); template_left(midpt-1:midpt+1, :) = 1;
template_right = template_left';

compare(template_left, template_right, temporal, 0, 'bars & decreasing');

%% Third test: sinusoid templates with flat temporal

os = linspace(0,pi,w);
[oo,~] = meshgrid(os,os);

temporal = linspace(1,0,frames)';
template_left = sin(oo);
template_right = template_left';

compare(template_left, template_right, temporal, 0, 'sinusoid & decreasing');



end