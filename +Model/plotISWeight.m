function plotISWeight(params)
%Model.PLOTISWEIGHT Plot importance-sampling weight as a function of the
%sampled x value and different running log-posterior values.

sig_x = sqrt(params.var_x);
p_match = params.p_match;

p_x_C = @(C, pC) [+C sig_x p_match*pC; -C sig_x (1-p_match)*pC];
p_x_Cp = p_x_C(+1, 1); % 'p' for plus
p_x_Cm = p_x_C(-1, 1); % 'm' for minus

xs = linspace(-2, 2, 1001);
log_ps = linspace(-1, 1, 5); % 10^log_p : 1 odds in favor of C=+1
ps = 1 ./ (exp(-log_ps) + 1);

is_weights = zeros(length(log_ps), length(xs));
update_Cp = mog.pdf(xs, p_x_Cp)';
update_Cm = mog.pdf(xs, p_x_Cm)';

for i=1:length(log_ps)
    post_C = ps(i);
    is_weights(i, :) = 1 ./ (mog.pdf(xs, p_x_Cp) * post_C + mog.pdf(xs, p_x_Cm) * (1 - post_C));
end

% norm = 1 ./ (mog.pdf(0, p_x_Cp)*0.5 + mog.pdf(0, p_x_Cm)*0.5);
% is_weights = is_weights/ norm;

figure;
subplot(2, 1, 1);
hold on;
leg = cell(1, length(log_ps));
for i=1:length(log_ps)
    plot(xs, is_weights(i, :));
    leg{i} = ['P(C+) = ' num2str(ps(i))];
end
xlabel('sampled x');
ylabel('IS weight');
legend(leg);
title('Importance-Sampling weight for different posteriors');

subplot(2, 1, 2);
hold on;
colors = lines(length(log_ps));
for i=1:length(log_ps)
    handles(i) = plot(xs, update_Cp .* is_weights(i, :), '-', 'Color', colors(i, :));
    plot(xs, update_Cm .* is_weights(i, :), '--', 'Color', colors(i, :));
end
xlabel('sampled x');
ylabel('P(x_i|C) * w_i');
legend(handles, leg);
title('Importance-Sampling full update');

end