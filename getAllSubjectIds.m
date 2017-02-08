function ids = getAllSubjectIds(datadir, expt_prefix)

id_file = fullfile(datadir, [expt_prefix 'subjectIDs.mat']);
contents = load(id_file);

ids = cellfun(@(subject) subject{2}, contents.id_map, 'UniformOutput', false);

end