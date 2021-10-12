function optim_results = fitNestedITBModels(base_params, signals, choices)
for iP=1:length(base_params)
    base_params.model = 'itb';
    base_params.updates = 1;
end

nPara = length(base_params);
nTrials = length(choices);
if iscell(choices)
    nTrials = sum(cellfun(@length, choices));
end

ibs_evals = round(sqrt(nTrials));

    function nll = negloglike(x, fields, distribs)
        [~, ll] = Fitting.choiceModelLogProb(Fitting.setParamsFields(base_params, fields, x), ...
            distribs, signals, choices, 1);
        nll = -ll;
    end

%% Round 1: fit unbounded, noiseless integrator
fields = {'prior_C', 'log_lapse', 'log_signal_scale'};
distribs = Fitting.defaultDistributions(fields, true);

LB = cellfun(@(f) distribs.(f).lb, fields);
UB = cellfun(@(f) distribs.(f).ub, fields);
psopts = optimoptions('patternsearch', 'display', 'final');
opts = optimoptions('fmincon', 'display', 'final');

for iInit=10:-1:1
    vals = cellfun(@(f) distribs.(f).priorrnd(1), fields);
    ps_min(iInit,:) = patternsearch(@(x) negloglike(x, fields, distribs), vals, ...
        [], [], [], [], LB, UB, [], psopts);
    [bestfit_vals(iInit,:), bestfit_nll(iInit)] = fmincon(@(x) negloglike(x, fields, distribs), ps_min(iInit,:), ...
        [], [], [], [], LB, UB, [], opts);
end

[~, iBest] = min(bestfit_nll);
optim_results{1} = Fitting.setParamsFields(base_params, fields, bestfit_vals(iBest, :));
[lp, ll, lv] = Fitting.choiceModelLogProb(optim_results{1}, distribs, signals, choices, [], 1);
for iP=1:nPara
    optim_results{1}(iP).fit_fields = fields;
    optim_results{1}(iP).fit_method = 'nested-ideal';
    optim_results{1}(iP).ll = ll;
end

%% Round 2: add a bound
fields = [fields {'log_bound'}];
distribs = Fitting.defaultDistributions(fields, true);
vals = Fitting.getParamsFields(optim_results{1}, fields);

slice = linspace(distribs.log_bound.plb, distribs.log_bound.pub);
for iSlc=length(slice):-1:1
    this_vals = vals;
    this_vals(strcmpi(fields, 'log_bound')) = slice(iSlc);
    slice_nll(iSlc) = negloglike(this_vals, fields, distribs);
end
[~,iMaxLL] = min(slice_nll);
vals(strcmpi(fields, 'log_bound')) = slice(iMaxLL);

LB = cellfun(@(f) distribs.(f).lb, fields);
UB = cellfun(@(f) distribs.(f).ub, fields);
bestfit_vals = fmincon(@(x) negloglike(x, fields, distribs), ...
    vals, [], [], [], [], LB, UB, [], opts);

optim_results{2} = Fitting.setParamsFields(base_params, fields, bestfit_vals);
[lp, ll, lv] = Fitting.choiceModelLogProb(optim_results{2}, distribs, signals, choices, [], 1);
for iP=1:nPara
    optim_results{2}(iP).fit_fields = fields;
    optim_results{2}(iP).fit_method = 'nested-itb';
    optim_results{1}(iP).ll = ll;
end

%% Round 3a: add a leak ('gamma')
fields = [fields {'gamma'}];
distribs = Fitting.defaultDistributions(fields, true);
vals = Fitting.getParamsFields(optim_results{2}, fields);

slice = linspace(distribs.gamma.plb, distribs.gamma.pub);
for iSlc=length(slice):-1:1
    this_vals = vals;
    this_vals(strcmpi(fields, 'gamma')) = slice(iSlc);
    slice_nll(iSlc) = negloglike(this_vals, fields, distribs);
end
[~,iMaxLL] = min(slice_nll);
vals(strcmpi(fields, 'gamma')) = slice(iMaxLL);

LB = cellfun(@(f) distribs.(f).lb, fields);
UB = cellfun(@(f) distribs.(f).ub, fields);
bestfit_vals = fmincon(@(x) negloglike(x, fields, distribs), ...
    vals, [], [], [], [], LB, UB, [], opts);

optim_results{3} = Fitting.setParamsFields(base_params, fields, bestfit_vals);
[lp, ll, lv] = Fitting.choiceModelLogProb(optim_results{3}, distribs, signals, choices, [], 1);
for iP=1:nPara
    optim_results{3}.fit_fields = fields;
    optim_results{3}.fit_method = 'nested-itb-leak';
    optim_results{3}(iP).ll = ll;
end

%% Round 3b: repeat 'gamma' but allow it in [-1 1]
fields{strcmpi(fields, 'gamma')} = 'neggamma';
distribs = Fitting.defaultDistributions(fields, true);
vals = Fitting.getParamsFields(optim_results{2}, fields);

slice = linspace(distribs.neggamma.plb, distribs.neggamma.pub);
for iSlc=length(slice):-1:1
    this_vals = vals;
    this_vals(strcmpi(fields, 'neggamma')) = slice(iSlc);
    slice_nll(iSlc) = negloglike(this_vals, fields, distribs);
end
[~,iMaxLL] = min(slice_nll);
vals(strcmpi(fields, 'neggamma')) = slice(iMaxLL);

LB = cellfun(@(f) distribs.(f).lb, fields);
UB = cellfun(@(f) distribs.(f).ub, fields);
bestfit_vals = fmincon(@(x) negloglike(x, fields, distribs), ...
    vals, [], [], [], [], LB, UB, [], opts);

optim_results{4} = Fitting.setParamsFields(base_params, fields, bestfit_vals);
[lp, ll, lv] = Fitting.choiceModelLogProb(optim_results{4}, distribs, signals, choices, [], 1);
for iP=1:nPara
    optim_results{4}(iP).fit_fields = fields;
    optim_results{4}(iP).fit_method = 'nested-itb-gamma';
    optim_results{3}(iP).ll = ll;
end

%% For each of the above deterministic models, now add a noise term and do stochastic search
for iModel=1:length(optim_results)
    det_base_model = optim_results{iModel};
    fields = [det_base_model.fit_fields {'log_noise'}];
    distribs = Fitting.defaultDistributions(fields, true);
    vals = Fitting.getParamsFields(det_base_model, fields);
    
    % Round 1: take a step along each parameter's marginal while building a table of data to train a
    % full joint GP. Go in reverse order to start with a step on the noise term.
    slice_size = 50;
    eval_points = zeros(length(fields)*slice_size, length(fields));
    eval_lls = zeros(length(fields)*slice_size, 1);
    eval_vars = zeros(length(fields)*slice_size, 1);
    for iPara=length(fields):-1:1
        slice = linspace(distribs.(fields{iPara}).plb, distribs.(fields{iPara}).pub, slice_size)';
        this_vals = repmat(vals, slice_size, 1);
        this_vals(:, iPara) = slice;
        slice_ll = zeros(size(slice));
        slice_var = zeros(size(slice));
        for iSlc=1:length(slice)
            [~, slice_ll(iSlc), slice_var(iSlc)] = Fitting.choiceModelLogProbIBS(...
                Fitting.setParamsFields(base_params, fields, this_vals(iSlc, :)), ...
                distribs, signals, choices, [], ibs_evals);
        end
        
        % Store this 1D slice of evaluations
        eval_points((iPara-1)*slice_size+1:iPara*slice_size, :) = this_vals;
        eval_lls((iPara-1)*slice_size+1:iPara*slice_size) = slice_ll;
        eval_vars((iPara-1)*slice_size+1:iPara*slice_size) = slice_var;
        
        % For this first round, take a coordinate ascent step using a GP on this parameter's
        % marginal slice alone; no GP fitting is done here.
        gp_scale(iPara) = (distribs.(fields{iPara}).pub - distribs.(fields{iPara}).plb) / 5;
        gp_sigma(iPara) = sqrt(max(1, var(slice_ll) - mean(slice_var)));
        gp_marg_ll = fitrgp(slice, slice_ll, 'Sigma', sqrt(mean(slice_var)), 'FitMethod', 'none', ...
            'KernelFunction', 'SquaredExponential', 'KernelParameters', [gp_scale(iPara) gp_sigma(iPara)], ...
            'Basis', 'PureQuadratic', 'Beta', negquadfit(slice, slice_ll));
        [~, vals(iPara)] = gpcoordinatestep(gp_marg_ll, vals(iPara), 1, min(slice), max(slice));
        
        subplot(1,length(fields),iPara); cla; hold on;
        errorbar(slice, slice_ll, sqrt(slice_var), '.');
        beta = negquadfit(slice, slice_ll);
        plot(slice, [ones(size(slice)) slice slice.^2]*beta)
        [y,ysd] = gp_marg_ll.predict(slice);
        errorbar(slice, y, ysd, '.k');
        plot(vals(iPara)*[1 1], ylim, '--r');
        title(sprintf('%s [marginal]', fields{iPara}));
        drawnow;
        
        % Debug printout
        [pred_ll, pred_ll_sd] = gp_marg_ll.predict(vals(iPara));
        fprintf('Round 1 Marginal Slice [%s] update [%f] -> [%f]. LL = %.2f+/-%.2f\n', ...
            fields{iPara}, Fitting.getParamsFields(det_base_model, fields{iPara}), vals(iPara), pred_ll, pred_ll_sd);
    end

    % Round 2: construct a GP over the joint data from all slices so far. See
    % @Fitting.marginalJointKernel for details on kernel params.
    kernel_params = [gp_scale; gp_sigma];
    kernel_params = [kernel_params(:); sqrt(mean(gp_sigma.^2))];
    kernel_params = log(kernel_params);

    [~, best_ll, best_ll_var] = Fitting.choiceModelLogProbIBS(Fitting.setParamsFields(base_params, fields, vals), ...
        distribs, signals, choices, [], ibs_evals);
    best_ll_sd = sqrt(best_ll_var);
    
    % Continue taking coordinate ascent steps using joint GP to convergence.
    tolx = cellfun(@(f) (distribs.(f).pub - distribs.(f).plb)/1000, fields);
    beta = negquadfit(eval_points, eval_lls);
    for itr=2:100
        vals(itr, :) = vals(itr-1, :);
        for iPara=length(fields):-1:1
            slice = linspace(distribs.(fields{iPara}).plb, distribs.(fields{iPara}).pub, slice_size)';
            this_vals = repmat(vals(itr, :), slice_size, 1);
            this_vals(:, iPara) = slice;
            slice_ll = zeros(size(slice));
            slice_var = zeros(size(slice));
            for iSlc=1:length(slice)
                [~, slice_ll(iSlc), slice_var(iSlc)] = Fitting.choiceModelLogProbIBS(...
                    Fitting.setParamsFields(base_params, fields, this_vals(iSlc, :)), ...
                    distribs, signals, choices, [], ibs_evals);
            end
            
            % Update table
            eval_points(end+1:end+slice_size, :) = this_vals;
            eval_lls(end+1:end+slice_size) = slice_ll;
            eval_vars(end+1:end+slice_size) = slice_var;
            
            % Update Joint GP
            [beta, resid] = negquadfit(eval_points, eval_lls, beta);
            kernel_params(2*iPara) = log(sqrt(max(1, var(slice_ll)-mean(slice_var))));
            joint_var = max(1, var(resid) - sum(exp(kernel_params(2:2:end-1))));
            kernel_params(end) = log(sqrt(joint_var));
            gp_joint_ll = fitrgp(eval_points, eval_lls, 'Sigma', sqrt(mean(eval_vars)), 'FitMethod', 'none', ...
                'KernelFunction', @Fitting.marginalJointKernel, 'KernelParameters', kernel_params, ...
                'Basis', 'PureQuadratic', 'Beta', beta);
            
            % GP coordinate ascent step
            [~, vals(itr, iPara)] = gpcoordinatestep(gp_joint_ll, vals(itr, :), iPara, min(slice), max(slice));

            % Debug printout
            [pred_ll, pred_ll_sd] = gp_joint_ll.predict(vals(itr, :));
            fprintf('Joint Iter %02d\tSlice %s\tupdate [%f] -> [%f]. LL = %.2f+/-%.2f\n', ...
                itr, fields{iPara}, vals(itr-1,iPara), vals(itr,iPara), pred_ll, pred_ll_sd);

            % Debug figure
            subplot(1,length(fields),iPara); cla; hold on;
            errorbar(slice, slice_ll, sqrt(slice_var), '.');
            beta = negquadfit(slice, slice_ll);
            plot(slice, [ones(size(slice)) slice slice.^2]*beta)
            [y,ysd] = gp_joint_ll.predict(this_vals);
            errorbar(slice, y, ysd, '.k');
            plot(vals(itr, iPara)*[1 1], ylim, '--r');
            title(sprintf('%s [joint]', fields{iPara}));
            drawnow;
        end
        
        % Check for convergence both in terms of vals (x) and estimated LL (y)
        delta_x = abs(vals(itr,:) - vals(itr-1,:));
        [best_ll(itr), best_ll_sd(itr)] = gp_joint_ll.predict(vals(itr, :));
        delta_ll = abs(best_ll(itr) - best_ll(itr-1)) / best_ll_sd(itr);
        
        if all(delta_x < tolx) && delta_ll < 0.1
            break;
        end
    end
    
    % Store this best 'noisy' version of the current model.
    optim_results{iModel} = Fitting.setParamsFields(base_params, fields, vals(end, :));
    [lp, ll, lv] = Fitting.choiceModelLogProbIBS(optim_results{iModel}, distribs, signals, choices, [], 10*ibs_evals);
    for iP=1:nPara
        optim_results{iModel}(iP).fit_fields = fields;
        optim_results{iModel}(iP).fit_method = 'nested-gp-ascent';
        optim_results{iModel}(iP).ll = ll;
    end
end
end

%%

function [beta, residuals] = negquadfit(x, y, beta0)
H = [ones(size(x,1),1) x x.^2];
if nargin < 3
    beta0 = [zeros(1, size(x,2)+1), -ones(1, size(x,2))]';
end
lb = -inf(size(beta0));
ub = +inf(size(beta0));
ub(size(x,2)+2:end) = 0;
opts = optimoptions('fmincon', 'Display', 'none');
beta = fmincon(@(b) sum((H*b-y).^2), beta0, [], [], [], [], lb, ub, [], opts);
residuals = y-H*beta;
end

function [maxval, maxx] = gpcoordinatestep(gp, startpoint, dim, lb, ub)
vals = repmat(startpoint, 1000, 1);
vals(:, dim) = linspace(lb, ub, 1000);
y = gp.predict(vals);
[maxval, imax] = max(y);
maxx = vals(imax);
end