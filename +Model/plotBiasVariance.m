function plotBiasVariance(trials, frames, prior, likelihood, sampling_params, gammas)
%PLOTBIASVARIANCE Run sampling model with different gammas to show bias
%variance tradeoff.

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

for i=1:length(gammas)
    gamma = gammas(i);
    params_copy = sampling_params;
    params_copy.gamma = gamma;
    [data, params_copy.var_e, prefix] = Model.genDataWithPriorLikelihood(trials, frames, prior, likelihood);
    [results, data] = Model.loadOrRunSamplingModel(data, prefix, params_copy);
    [ideal_results] = Model.runIdealObserver(data, params_copy);
    
    % Get difference between sampling log posterior and ideal log posterior
    log_posterior_difference = results.walk - ideal_results.walk;
    [M, L, U] = meanci(log_posterior_difference);
    
    stringID = Model.getModelStringID(prefix, params_copy);
    
    fig = figure();
    boundedline(1:frames+1, M', [M-L; U-M]');
    title(['Log Posterior Odds difference from Ideal ' strrep(stringID, '_', ' ')]);
    saveas(fig, fullfile(savedir, ['biasVariance_' stringID '.fig']));
end
end

