function fig = visModelFitPsychometric(SubjectData, phase, VP, fields, params, n_inner, n_samples)

switch phase
    case 1
        uVar = unique(SubjectData.true_ratio);
        trialfn = @(ratio) SubjectData.true_ratio == ratio;
        xlab = 'ratio';
    case 2
        uVar = unique(SubjectData.noise .* sign(SubjectData.correct_answer-.5));
        trialfn = @(kappa) SubjectData.noise == abs(kappa) & ...
            (sign(SubjectData.correct_answer-.5)==sign(kappa) | kappa == 0);
        xlab = 'kappa';
end

for iVar=length(uVar):-1:1
    trials = trialfn(uVar(iVar));
    [data_pm(iVar), data_pmci(:,iVar)] = binofit(sum(SubjectData.choice(trials) == +1), sum(trials));
end

Xsamp = vbmc_rnd(VP, n_samples);
for iSamp=n_samples:-1:1
    disp(iSamp);
    [~, model_params, data_translated] = Fitting.subjectDataLogLikelihood(...
        Xsamp(iSamp,:), fields, params, SubjectData, n_inner);
    results = Model.runVectorized(model_params, data_translated);
    for iVar=length(uVar):-1:1
        trials = trialfn(uVar(iVar));
        mdl_pm(iSamp,iVar) = binofit(sum(results.choices(trials) == +1), sum(trials));
    end
end

fig = figure; hold on;
plot(uVar, mdl_pm', 'Color', [.6 .6 .6]);
errorbar(uVar, data_pm, data_pm-data_pmci(1,:), data_pmci(2,:)-data_pm, '-ob');
xlabel(xlab);
ylabel('% chose rightward');
set(gca, 'YTick', 0:.1:1, 'YTickLabel', 0:10:100);
grid on;

end