function params = setParamsFields(params, fields, values)
%FITTING.SETPARAMSFIELDS helper / wrapper to set fields of 'params' to the given 'values' with some
%extra processing for cases where params is a struct array, where fields are prefixed with 'log_',
%suffixed with a number like 'lapse1', or where fields have an index like gamma_1.
%
%The inverse operations are done in @Fitting.getParamsFields
%
% Example 1: basic usage
%     params = Model.newModelParams();
%     disp([params.gamma params.lapse]);
%     params = Fitting.setParamsFields(params, {'gamma', 'lapse'}, [0.1 0.05]);
%     disp([params.gamma params.lapse]);
%
% Example 2: log_ prefix
%     params = Model.newModelParams();
%     disp([params.temperature]);
%     params = Fitting.setParamsFields(params, 'log_temperature', 0);
%     disp([params.temperature]);
%
% Example 3: parameter array
%     params = [Model.newModelParams() Model.newModelParams()];
%     disp([params.temperature]);
%     disp([params.gamma]);
%     params = Fitting.setParamsFields(params, {'log_temperature', 'gamma_1', 'gamma_2'}, [0 .1 .2]);
%     disp([params.temperature]);
%     disp([params.gamma]);

if ~iscell(fields), fields = {fields}; end
for iF=1:length(fields)
    % If, e.g. field is 'log_temperature', then we set params.temperature to exp() of the value
    if startsWith(fields{iF}, 'log_')
        values(iF) = exp(values(iF));
        fields{iF} = strrep(fields{iF}, 'log_', '');
    end

    % If the field is, e.g. gamma_1 or gamma_2, that indicates that this parameter should be
    % applied to only one of the array of parameters, e.g. params(1).gamma = gamma_1 and
    % params(2).gamma = gamma_2.
    hasindex = ~isempty(regexpi(fields{iF}, '\w+\d*_\d+$'));
    if hasindex
        % Set the value on just the struct with the selected index
        idxpart = regexprep(fields{iF}, '\w+\d*_(\d+)', '$1');
        namepart = regexprep(fields{iF}, '(\w+\d*)_\d+', '$1');
        idx = str2double(idxpart);
        params(idx).(namepart) = values(iF);
    else
        % Set the value on all 'params' structs
        for iP=1:length(params)
            params(iP).(fields{iF}) = values(iF);
        end
    end
end
params = arrayfun(@Fitting.sanitize, params);
end