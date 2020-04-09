function values = getParamsFields(params, fields)
%FITTING.GETPARAMSFIELDS inverse of @Fitting.setParamsFields. See comments there.
if ~iscell(fields), fields = {fields}; end
for iF=length(fields):-1:1
    hasindex = ~isempty(regexpi(fields{iF}, '\w+\d*_\d+$'));
    if hasindex
        idxpart = regexprep(fields{iF}, '\w+\d*_(\d+)', '$1');
    	idx = str2double(idxpart);
        fields{iF} = regexprep(fields{iF}, '(\w+\d*)_\d+', '$1');
    else
    	idx = 1;
    end

    if startsWith(fields{iF}, 'log_')
        values(iF) = log(params(idx).(strrep(fields{iF}, 'log_', '')));
    else
        values(iF) = params(idx).(fields{iF});
    end
end
end