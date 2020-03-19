function values = getParamsFields(params, fields)
if ~iscell(fields), fields = {fields}; end
for iF=length(fields):-1:1
    if startsWith(fields{iF}, 'log_')
        values(iF) = log(params.(strrep(fields{iF}, 'log_', '')));
    else
        values(iF) = params.(fields{iF});
    end
end
end