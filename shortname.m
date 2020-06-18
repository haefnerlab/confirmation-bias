function name = shortname(name)

if contains(name, 'subject')
    idx = regexpi(name, 'subject');
    idstart = idx + length('subject');
    name = ['s' name(idstart:end)];
else
    name = upper(name);
end

end