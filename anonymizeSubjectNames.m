function anonymizeSubjectNames

datadir = fullfile(pwd, '..', 'RawData');
datafiles = dir(fullfile(datadir, '*-subjectIDs.mat'));

for i=1:length(datafiles)
    file = datafiles(i);
    fullname = fullfile(datadir, file.name);
    load(fullname);
    disp(file.name);
    for j=1:length(id_map)
        name = id_map{j}{1};
        if isnumeric(id_map{j}{1})
            id_map{j}{1} = sprintf('%X', id_map{j}{1});
        elseif regexp(name, '[a-z]+')
            id_map{j}{1} = sprintf('%X', string2hash(name));
            fprintf(['\t' name '\t' id_map{j}{1} '\n']);
        end
    end
    save(fullname, 'id_map');
end