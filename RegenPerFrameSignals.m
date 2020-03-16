function perFrameSigs = RegenPerFrameSignals(SubjectData, kernelKappa)
%REGENPERFRAMESIGNALS Regenerate stimuli from seed and compute amount of signal as difference of
%signal from fourier-domain templates between 2 categories.

parfor tr=1:SubjectData.current_trial
    % Regenerate stimulus from seed (can be slow if error-recovery is needed)
    im = GaborStimulus(SubjectData, tr, 2);

    % Center the image after converting out of uint8 type
    im = double(im)-127;
    
    % Determine template kappa - if given kernelKappa is [], use the 'true' value for this trial,
    % minimum of .04. Else use the same template for all trials.
    if isempty(kernelKappa)
        thisKernelKappa = max(.04, SubjectData.noise(tr));
    else
        thisKernelKappa = kernelKappa;
    end
    perFrameSigs(tr, :) = bpg.getSignal(im, SubjectData.left_category, thisKernelKappa) - ...
        bpg.getSignal(im, SubjectData.right_category, thisKernelKappa);
end

end