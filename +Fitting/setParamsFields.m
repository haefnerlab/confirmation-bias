function params = setParamsFields(params, fields, values)
if ~iscell(fields), fields = {fields}; end
for iF=1:length(fields)
    if startsWith(fields{iF}, 'log_')
        values(iF) = exp(values(iF));
        fields{iF} = strrep(fields{iF}, 'log_', '');
    end
	for iP=1:length(params)
	    params(iP).(fields{iF}) = values(iF);
	end
end
params = Fitting.sanitize(params);
end