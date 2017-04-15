function fig=Get_Figure(tag,option,idx,prop)

if nargin<2, option='clear'; end
if nargin<3, idx=1; end
if nargin<4, prop=''; end

if isempty(findstr(option,'intact'))
    option=[option 'clear'];
end

figs=findobj('Tag',tag);

if isempty(figs) | (~isempty(figs) & length(figs)<idx) | ~isempty(findstr(option,'force'))
    switch length(prop)
        case 0, fig=figure;
        case 2, fig=figure(prop{1},prop{2});
    end
else
    fig=figs(idx);
    figure(fig(idx));
end

if ~isempty(findstr(option,'clear'))
    clf(fig);
end

set(fig,'Tag',tag);
hold on;

if ~(isempty(tag) | tag(1)==' ' | tag=='~')
  title(tag);
end
