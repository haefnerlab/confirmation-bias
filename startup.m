if exist('psignifit', 'file')~=2, addpath('psignifit-master'); end
if exist('smoothn', 'file')~=2, addpath('smoothn'); end
if exist('fadecolors', 'file')~=2, addpath('tools'); end
if exist('boundedline', 'file')~=2
    D = dir(fullfile('tools', 'boundedline-pkg'));
    for i=1:length(D)
        if D(i).isdir
            addpath(fullfile('tools', 'boundedline-pkg', D(i).name));
        end
    end
end
