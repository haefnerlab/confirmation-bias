function anonId = getSubjectId(data_directory, prefix)
if nargin < 2, prefix = ''; end

prompts = {'Enter subject name:', '[N]ew or [R]eturning?'};
answers = inputdlg(prompts);
[subject_name, new_or_return] = answers{:};
subject_name = lower(subject_name);

if ~any(strcmpi(new_or_return, {'n', 'r'}))
    error('Enter ''n'' or ''r'' for %s', prompts{2});
end

data_file = fullfile(data_directory, [prefix 'subjectIDs.mat']);

if ~exist(data_file, 'file')
    id_map = {};
else
    load(data_file);
end

for i=1:length(id_map)
    if strcmpi(id_map{i}{1}, subject_name)
        if strcmpi(new_or_return, 'n')
            error('%s is returning?', subject_name);
        end
        anonId = id_map{i}{2};
        return
    end
end

% If reached here, subject_name is not in our database.
if strcmpi(new_or_return, 'r')
    error('%s is new?', subject_name);
end

% Create new entry.
anonId = sprintf('%ssubject%02d', prefix, length(id_map)+1);
id_map{end+1} = {subject_name, anonId};
save(data_file, 'id_map');
end