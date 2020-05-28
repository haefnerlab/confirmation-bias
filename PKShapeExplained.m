function [partialBetas, partialBetaErrs, ablatedFields] = PKShapeExplained(signals, params, test_fields, all_fields, values)

assert(length(params) == 1, 'Only makes sense on / only supports 1 condition at a time');

if nargin < 5, values = Fitting.getParamsFields(params, all_fields); end

% Per row of 'values' (e.g. if sampled), estimate PK shape on the given data with each of the given
% 'test_fields' and each combination of them returned to defaults. This is 2^|test_fields|
% evaluations per row
ablatedFields = {{}};
for iF=1:length(test_fields)
    extendedFields = cellfun(@(otherFields) horzcat(otherFields, test_fields{iF}), ablatedFields, 'uniformoutput', false);
    ablatedFields = horzcat(ablatedFields, extendedFields);
end

partialBetas = zeros(size(values,1), length(ablatedFields));
partialBetaErrs = zeros(size(values,1), length(ablatedFields));
for iVal=1:size(values, 1)
    for iAbl=1:length(ablatedFields)
        % Begin this run by setting all fields to values(iVal, :)
        this_params = Fitting.setParamsFields(params, all_fields, values(iVal, :));

        % Now return all fields to be ablated back to their null/default state
        for iF=1:length(ablatedFields{iAbl})
            if contains(ablatedFields{iAbl}{iF}, 'gamma')
                this_params.gamma = 0;
            elseif contains(ablatedFields{iAbl}{iF}, 'bound')
                this_params.bound = inf;
            elseif contains(ablatedFields{iAbl}{iF}, 'samples')
                this_params.samples = 100;
            elseif contains(ablatedFields{iAbl}{iF}, 'noise')
                this_params.noise = 0;
            else
                error('Unsupported field: %s', ablatedFields{iAbl}{iF});
            end
        end

        % Using these parameters, simulate choices and estimate PK shape
        results = Model.runVectorized(this_params, signals/this_params.signal_scale);
        [abb, ~, errors] = CustomRegression.ExponentialPK(signals/this_params.signal_scale, results.choices == +1, false);
        partialBetas(iVal, iAbl, :) = abb(2);
        partialBetaErrs(iVal, iAbl, :) = errors(2);
    end
end

end