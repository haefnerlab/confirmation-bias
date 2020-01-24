function Data = GaborThresholdTrials(Data, phase, threshold, floor)

if phase == 0
    stair_param = 'contrast';
elseif phase == 1
    stair_param = 'true_ratio';
elseif phase == 2
    stair_param = 'noise';
end

test_trials = Data.(stair_param) <= threshold;
if exist('floor', 'var')
    test_trials = test_trials & Data.(stair_param) >= floor;
end

%% Subselect trial data.

Data.current_trial = sum(test_trials);
Data = maskFieldsDetectDimension(Data, test_trials);

end

function str = maskFieldsDetectDimension(str, mask)
% Apply logical mask to fields in the given struct, automatically detecting which fields and which
% dimension of those fields need processing based on their length.
nTrials = length(mask);
fields = fieldnames(str);
for iField=1:length(fields)
    needsMask = false;
    val = str.(fields{iField});
    dims = ndims(val);
    ref = struct('type', '()', 'subs', {repmat({':'}, 1, dims)});
    for d=1:dims
        if size(val, d) == nTrials
            ref.subs{d} = mask;
            needsMask = true;
            break
        end
    end
    if needsMask
        str.(fields{iField}) = subsref(val, ref);
    end
end
end
