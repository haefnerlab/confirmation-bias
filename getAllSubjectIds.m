function ids = getAllSubjectIds(datadir, expt_prefix, phase)

id_file = fullfile(datadir, [expt_prefix 'subjectIDs.mat']);
contents = load(id_file);

ids = cellfun(@(subject) subject{2}, contents.id_map, 'UniformOutput', false);

if exist('phase', 'var')
    switch phase
        case 0, exptType = 'Contrast';
        case 1, exptType = 'Ratio';
        case 2, exptType = 'Noise';
        otherwise, error('phase must be 0, 1, or 2'); 
    end
    mask = true(size(ids));
    for i=1:length(ids)
        files = dir(fullfile(datadir, [ids{i} '*' exptType '.mat']));
        filesQuit = dir(fullfile(datadir, [ids{i} '*' exptType 'Quit.mat']));
        if isempty(files) && isempty(filesQuit)
            mask(i) = false;
        end
    end
    ids = ids(mask);
end

end