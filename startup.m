if exist('psignifit', 'file')~=2, addpath('psignifit-master'); end
if exist('smoothn', 'file')~=2, addpath('smoothn'); end
if exist('fadecolors', 'file')~=2, addpath('tools'); end
if exist('advancedcolormap', 'file')~=2, addpath('AdvancedColorMap'); end
if exist('vbmc', 'file')~=2, addpath(genpath('vbmc')); end
if exist('bads', 'file')~=2, addpath(genpath('bads')); end
if exist('boundedline', 'file')~=2
    addpath(genpath('tools/boundedline-pkg'));
end
if exist('processManager', 'file')~=2, addpath('procman'); end

%% Create useful workspace variables for analysis

RATIO_PHASE = 1; % aka HSLC
NOISE_PHASE = 2; % aka LSHC
THRESHOLD = 0.7;
KERNEL_KAPPA = 0.16;
DATADIR = fullfile(pwd, '..', 'PublishData');
MEMODIR = fullfile(pwd, '..', 'Precomputed');

DARK_RED = [204 0 0] / 255;
DARK_BLUE = [32 74 135] / 255;
% fade 50% towards white
LIGHT_RED = DARK_RED * .5 + .5;
LIGHT_BLUE = DARK_BLUE * .5 + .5;

ratioSubjects = arrayfun(@(i) sprintf('BPGTask-subject%02d', i), setdiff(1:15, [1 5 12]), 'UniformOutput', false);
noiseSubjects = arrayfun(@(i) sprintf('BPGTask-subject%02d', i), setdiff(1:15, [1 5]), 'UniformOutput', false);

% Informed is the opposite of naive. Subjects 7 through 9 were authors.
informedSubjects = arrayfun(@(i) sprintf('BPGTask-subject%02d', i), [7 8 9], 'UniformOutput', false);

naiveRatioSubjects = setdiff(ratioSubjects, informedSubjects);
naiveNoiseSubjects = setdiff(noiseSubjects, informedSubjects);
naiveBothSubjects = intersect(naiveRatioSubjects, naiveNoiseSubjects);

% Re-order so that informed subjects are at the end
ratioSubjects = [naiveRatioSubjects informedSubjects];
noiseSubjects = [naiveNoiseSubjects informedSubjects];
bothSubjects = sort([intersect(naiveRatioSubjects, naiveNoiseSubjects) informedSubjects]);
is_naive = ismember(bothSubjects, naiveBothSubjects);

% Compute population-level psycho metrics
N = length(bothSubjects);

for n=N:-1:1
    hyphenSplit = strsplit(bothSubjects{n}, '-');
    shortnames{n} = hyphenSplit{2};
end