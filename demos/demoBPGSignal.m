function demoBPGSignal(center)

if nargin < 1, center = 0; end

oris = linspace(0,180);

GaborData = newGaborData; % Use defaults

signals = [22.5 45 90];
leg = arrayfun(@num2str, signals, 'UniformOutput', false);

figure; hold on;

for s=signals
[im, ~] = bpg.genImages(length(oris), GaborData.stim_size, GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, oris, s, GaborData.annulus);
emp_signals = bpg.getSignal(im, center, s);
plot(oris, emp_signals);
end

xlabel('rotation from true orientation');
ylabel('fourier signal');
title('Signal level vs offset');
legend(leg);

end