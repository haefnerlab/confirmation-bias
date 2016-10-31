function isFieldResult = myIsField(inStruct, fieldName)
% A function to check for the existence of a field in a struct
% inStruct is the name of the structure or an array of structures to search
% fieldName is the string name of the field for which the function searches

isFieldResult = 0;
f = fieldnames(inStruct(1));
for i = 1:length(f)
    if(strcmp(f{i}, strtrim(fieldName)))
        isFieldResult = 1;
        return;
    elseif isstruct(inStruct(1).(f{i}))
        isFieldResult = myIsField(inStruct(1).(f{i}), fieldName);
        if isFieldResult
            return;
        end
    end
end