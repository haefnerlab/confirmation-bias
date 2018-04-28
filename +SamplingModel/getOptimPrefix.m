function prefix = getOptimPrefix(variables, grid_search_size)
if isempty(variables)
    prefix = '';
    return
end

prefix = 'optim';

if nargin > 1 && ~isempty(grid_search_size)
    prefix = [prefix '_grid' num2str(grid_search_size)];
end

if any(strcmpi('p_match', variables))
    prefix = [prefix '_PM'];
end
if any(strcmpi('var_s', variables))
    prefix = [prefix '_VE'];
end
if any(strcmpi('gamma', variables))
    prefix = [prefix '_G'];
end
if any(strcmpi('prior_C', variables))
    prefix = [prefix '_PC'];
end
end

