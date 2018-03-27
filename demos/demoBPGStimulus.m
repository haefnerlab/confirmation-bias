function demoBPGStimulus

sz = 100;
oriDEG = 0;
oriKappa = 0.3;
spFreqCycles = 6;
spFreqStdCycles = 1;
annulus = 10;
seed = randseed;

f = figure;

    function genAndPlot()
        rng(seed, 'twister');
        [im, imF] = bpg.genImages(1, sz, spFreqCycles / sz, spFreqStdCycles / sz, ...
            oriDEG, oriKappa, annulus);
        emp_imF = fftshift(fft2(squeeze(im)));
        ax1 = subplot(1, 3, 1);
        imagesc(squeeze(im));
        colormap gray;
        axis image;
        set(ax1, 'XTick', []);
        set(ax1, 'YTick', []);
        set(ax1, 'YDir', 'reverse');
        colorbar;
        title('Stimulus');

        ax2 = subplot(1, 3, 2);
        imagesc(abs(squeeze(imF)));
        colormap gray;
        axis image;
        set(ax2, 'XTick', []);
        set(ax2, 'YTick', []);
        set(ax2, 'YDir', 'reverse');
        colorbar;
        title('Fourier Domain');

        ax2 = subplot(1, 3, 3);
        imagesc(abs(squeeze(emp_imF)));
        colormap gray;
        axis image;
        set(ax2, 'XTick', []);
        set(ax2, 'YTick', []);
        set(ax2, 'YDir', 'reverse');
        colorbar;
        title('Fourier Domain (with aperture)');
    end

    function randomizeSeed(varargin)
    seed = randseed;
    genAndPlot();
    end

genAndPlot();

%% Create ui controls
figPos = f.Position;

sliders(1) = uicontrol('Parent', f, 'Style', 'slider', 'min', 0, 'max', 180, 'value', oriDEG, ...
    'Callback', @(es, ed) sliderUpdate('oriDEG', es.Value));
labels{1} = 'Orientation';

sliders(2) = uicontrol('Parent', f, 'Style', 'slider', 'min', 0, 'max', 2.5, 'value', oriKappa, ...
    'Callback', @(es, ed) sliderUpdate('oriKappa', es.Value));
labels{2} = 'Ori. Kappa';

sliders(3) = uicontrol('Parent', f, 'Style', 'slider', 'min', 0, 'max', 20, 'value', spFreqCycles, ...
    'Callback', @(es, ed) sliderUpdate('spFreqCycles', es.Value));
labels{3} = 'Sp. Freq';

sliders(4) = uicontrol('Parent', f, 'Style', 'slider', 'min', 0, 'max', 20, 'value', spFreqStdCycles, ...
    'Callback', @(es, ed) sliderUpdate('spFreqStdCycles', es.Value));
labels{4} = 'Sp. F. Width';

sliders(5) = uicontrol('Parent', f, 'Style', 'slider', 'min', 0, 'max', 20, 'value', annulus, ...
    'Callback', @(es, ed) sliderUpdate('annulus', es.Value));
labels{5} = 'annulus';

buttonsize = [50, 20];
uicontrol('Parent', f, 'Style', 'pushbutton', 'String', 'New Seed', ...
    'Callback', @randomizeSeed, 'Position', [figPos(3)/2-buttonsize(1)/2, 65, buttonsize(1), buttonsize(2)]);

sliderWidth = figPos(3) / length(sliders);
for i=1:length(sliders)
    sliders(i).Position = [sliderWidth*(i-1+.05), 20, sliderWidth*.9, 20];
    uicontrol('Parent', f, 'Style', 'text', 'String', ...
        [num2str(sliders(i).Min) ' <-- ' labels{i} ' --> ' num2str(sliders(i).Max)], ...
        'Position', [sliderWidth*(i-1+.05), 40, sliderWidth*.9, 20]);
end

%% Updater function
    function sliderUpdate(param, value)
    eval([param ' = ' num2str(value) ';']);
    genAndPlot();
    end
end