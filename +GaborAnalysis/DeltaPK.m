function [perSubjectFigs, combinedfig] = DeltaPK(subjectIDs, phases, datadir)
%GABORANALYSIS.DELTAPK creates one figure per subject and a combined figure
%(if 2 or more subjects) showing temporal psychophysical kernel analysis.
%
% GABORANALYSIS.DELTAPK(subjectIDs, phases, datadir)
%
% Inputs:
% - subjectIDs: a cell array of subject IDs (strings)
% - phases: 0 for contrast, 1 for ratio, 2 for noise. If 2 numbers given,
%           the difference phases(1)-phases(2) is plotted.
% - datadir: (optional) override the default place to look for data files.

if nargin < 4, datadir = fullfile(pwd, '..', 'RawData'); end

catdir = fullfile(datadir, '..', 'ConcatData');
if ~exist(catdir, 'dir'), mkdir(catdir); end

memodir = fullfile(datadir, '..', 'Precomputed');
if ~exist(memodir, 'dir'), mkdir(memodir); end

    function [M, L, U, trials, frames, SubjectDataThresh] = getSubjectKernel(subjectId, phase)
        stair_var = get_stair_var(phase);
        SubjectData = LoadOrRun(@LoadAllSubjectData, ...
            {subjectId, phase, datadir}, fullfile(catdir, [subjectId '-' stair_var '.mat']));
        [floor, thresh] = GaborAnalysis.getThresholdWindow(subjectId, phase, .6, .75, datadir);
        trials = SubjectData.(stair_var) <= thresh & SubjectData.(stair_var) >= floor;
        SubjectDataThresh = GaborThresholdTrials(SubjectData, phase, thresh, floor);
        memo_name = ['Boot-PK-' stair_var '-' subjectId '-' num2str(thresh) '-' num2str(floor) '.mat'];
        [M, L, U] = LoadOrRun(@BootstrapWeightsGabor, ...
            {SubjectDataThresh, 500}, ...
            fullfile(memodir, memo_name));
        frames = SubjectData.number_of_images;
    end

CombinedSubjectDataByPhase = cell(size(phases));
CombinedIdentifier = '';

for i=1:length(subjectIDs)
    perSubjectFigs(i) = figure;
    subjectId = subjectIDs{i};
    if i == 1
        CombinedIdentifier = subjectId;
    else
        CombinedIdentifier = [CombinedIdentifier '-' subjectId];
    end
    if length(phases) == 1
        [M, L, U, trials, frames, SubjectDataThresh] = getSubjectKernel(subjectId, phases);
        boundedline(1:frames, M(1:frames)', [U(1:frames)-M(1:frames); M(1:frames)-L(1:frames)]');
        title([strrep(get_stair_var(phases), '_', ' ') 'temporal kernel (' num2str(sum(trials)) '/' num2str(length(trials)) ')']);
        xlim([-inf, inf]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', 0);
        xlabel('Time');
        ylabel('Psychophysical Kernel');
        
        if i == 1
            CombinedSubjectDataByPhase{1} = SubjectDataThresh;
        else
            CombinedSubjectDataByPhase{1} = ConcatGaborData(CombinedSubjectDataByPhase{1}, SubjectDataThresh);
        end
    else
        % 2 subplots: 2 pks in one and their difference in the other
        subplot(1, 2, 1);
        hold on;
        [M1, L1, U1, trials1, frames1, SubjectDataThresh1] = getSubjectKernel(subjectId, phases(1));
        %         leg{1} = [strrep(get_stair_var(phases(1)), '_', ' ') ' (' num2str(sum(trials1)) '/' num2str(length(trials1)) ')'];
        leg{1} = ['Task A (' num2str(sum(trials1)) '/' num2str(length(trials1)) ')'];
        
        [M2, L2, U2, trials2, frames2, SubjectDataThresh2] = getSubjectKernel(subjectId, phases(2));
        %         leg{2} = [strrep(get_stair_var(phases(2)), '_', ' ') ' (' num2str(sum(trials2)) '/' num2str(length(trials2)) ')'];
        leg{2} = ['Task B (' num2str(sum(trials2)) '/' num2str(length(trials2)) ')'];
        
        h = boundedline(1:frames1, M1(1:frames1)', [U1(1:frames1)-M1(1:frames1); M1(1:frames1)-L1(1:frames1)]', 'r', ...
            1:frames2, M2(1:frames2)', [U2(1:frames2)-M2(1:frames2); M2(1:frames2)-L2(1:frames2)]', 'b', ...
            'alpha');
        title('temporal kernels');
        xlim([-inf, inf]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', 0);
        xlabel('Time');
        ylabel('Psychophysical Kernel');
        legend(h, leg, 'Location', 'best');
        
        subplot(1, 2, 2);
        assert(frames1 == frames2, 'Cannot subtract PKs with different # frames');
        Md = M1 - M2;
        Ld = L1 - U2;
        Ud = U1 - L2;
        boundedline(1:frames1, Md(1:frames1)', [Ud(1:frames1)-Md(1:frames1); Md(1:frames1)-Ld(1:frames1)]');
        title('temporal kernel difference');
        xlim([-inf, inf]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', 0);
        xlabel('Time');
        ylabel('Psychophysical Kernel');
        
        if i == 1
            CombinedSubjectDataByPhase{1} = SubjectDataThresh1;
            CombinedSubjectDataByPhase{2} = SubjectDataThresh2;
        else
            CombinedSubjectDataByPhase{1} = ConcatGaborData(CombinedSubjectDataByPhase{1}, SubjectDataThresh1);
            CombinedSubjectDataByPhase{2} = ConcatGaborData(CombinedSubjectDataByPhase{2}, SubjectDataThresh2);
        end
    end
    
    perSubjectFigs(i).PaperUnits = 'inches';
    perSubjectFigs(i).PaperSize = [8 4];
    perSubjectFigs(i).PaperPosition = [0 0 8 4];
    saveas(perSubjectFigs(i), [subjectId '-PKPlot.png']);
end

combinedfig = -1;
if length(subjectIDs) > 1
    combinedfig = figure;
    
    if length(phases) == 1
        [M, L, U] = BootstrapWeightsGabor(CombinedSubjectDataByPhase{1}, 500);
        frames = CombinedSubjectDataByPhase{1}.number_of_images;
        boundedline(1:frames, M(1:frames)', [U(1:frames)-M(1:frames); M(1:frames)-L(1:frames)]');
        title([strrep(get_stair_var(phases), '_', ' ') ' combined temporal kernel']);
        xlim([-inf, inf]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', 0);
        xlabel('Time');
        ylabel('Psychophysical Kernel');
    else
        % 2 subplots: 2 pks in one and their difference in the other
        subplot(1, 2, 1);
        hold on;
        [M1, L1, U1] = BootstrapWeightsGabor(CombinedSubjectDataByPhase{1}, 500);
        frames1 = CombinedSubjectDataByPhase{1}.number_of_images;
        %         leg{1} = strrep(get_stair_var(phases(1)), '_', ' ');
        leg{1} = 'Task A';
        
        [M2, L2, U2] = BootstrapWeightsGabor(CombinedSubjectDataByPhase{2}, 500);
        frames2 = CombinedSubjectDataByPhase{2}.number_of_images;
        %         leg{2} = strrep(get_stair_var(phases(2)), '_', ' ');
        leg{2} = 'Task B';
        
        h = boundedline(1:frames1, M1(1:frames1)', [U1(1:frames1)-M1(1:frames1); M1(1:frames1)-L1(1:frames1)]', 'r', ...
            1:frames2, M2(1:frames2)', [U2(1:frames2)-M2(1:frames2); M2(1:frames2)-L2(1:frames2)]', 'b', ...
            'alpha');
        title('combined temporal kernels');
        xlim([-inf, inf]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', 0);
        xlabel('Time');
        ylabel('Psychophysical Kernel');
        legend(h, leg, 'Location', 'best');
        
        subplot(1, 2, 2);
        assert(frames1 == frames2, 'Cannot subtract PKs with different # frames');
        Md = M1 - M2;
        Ld = L1 - U2;
        Ud = U1 - L2;
        boundedline(1:frames1, Md(1:frames1)', [Ud(1:frames1)-Md(1:frames1); Md(1:frames1)-Ld(1:frames1)]');
        title('combined temporal kernel difference');
        xlim([-inf, inf]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', 0);
        xlabel('Time');
        ylabel('Psychophysical Kernel');
    end
    
    combinedfig.PaperUnits = 'inches';
    combinedfig.PaperSize = [8 4];
    combinedfig.PaperPosition = [0 0 8 4];
    saveas(combinedfig, [CombinedIdentifier '-PKPlot.png']);
end
end

function stair_var = get_stair_var(phase)
if phase == 0
    stair_var = 'contrast';
elseif phase == 1
    stair_var = 'true_ratio';
elseif phase == 2
    stair_var = 'noise';
else
    error('Expected phase 0 for Contrast or 1 for Ratio or 2 for Noise');
end
end