function [weights, errors, fig] = plotPK(params, pk_hprs)
%VARIATIONALMODEL.PLOTPK(params) plot PK of VB model for given params.
%
% [weights, errors, fig] = VARIATIONALMODEL.PLOTPK(params) returns PK and fig handle

savedir = fullfile('+VariationalModel', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

if nargin < 2, pk_hprs = [1 0 10]; end

results_uid = VariationalModel.getModelStringID(params);    
[~, data, choices] = LoadOrRun(@VariationalModel.runLatentZNormalX, {params}, ...
    fullfile(params.save_dir, results_uid));

% Randomly flip trial signs
[data, choices] = flipTrials(data, choices);

% Do regression
[weights, ~, errors] = CustomRegression.PsychophysicalKernel(data, choices, ...
    pk_hprs(1), pk_hprs(2), pk_hprs(3));

pk_id = ['PK_' results_uid];
savefile = fullfile(savedir, [pk_id '.fig']);

fig = figure(); hold on;
errorbar(1:params.frames, weights(1:end-1), errors(1:end-1));
errorbar(params.frames+1, weights(end), errors(end), '-r');
xlabel('time');
ylabel('weight');
saveas(fig, savefile);

end

function [data, choices] = flipTrials(data, choices)
flip_indexes = rand(length(choices), 1) < 0.5;
data(flip_indexes, :) = -data(flip_indexes, :);
choices = choices == +1;
choices(flip_indexes) = ~choices(flip_indexes);
end