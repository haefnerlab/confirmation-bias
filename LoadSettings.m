function settings = LoadSettings(path)

if nargin < 1, path = pwd; end
if isdir(path), path = fullfile(path, 'settings.mat'); end
if exist(path, 'file')
    disp(['Loading ', path]);
    contents = load(path);
    settings = contents.settings;
else
    disp('Using default settings');
    settings = struct();
end

% gammaTableFile should be a mat file containing a gamma table in a field
% whose name is in settings.gammaTable.
if ~isfield(settings, 'gammaTableFile'), settings.gammaTableFile = ''; end
if ~isfield(settings, 'gammaTable'), settings.gammaTable = ''; end
% PTB settings
if ~isfield(settings, 'whichScreen'), settings.whichScreen = 0; end
if ~isfield(settings, 'useOpenGL'), settings.useOpenGL = true; end
if ~isfield(settings, 'ptbSkipSyncTests'), settings.ptbSkipSyncTests = 0; end
trueSize = Screen('Resolution', settings.whichScreen);
if ~isfield(settings, 'screenSize'), settings.screenSize = [0 0 trueSize.width trueSize.height]; end
settings.monitorFPS = Screen('NominalFrameRate', settings.whichScreen, 1);
if settings.monitorFPS == 0
    warning('Getting frame rate failed. Assuming exactly 60.0 FPS');
    settings.monitorFPS = 60.0;
end
if ~isfield(settings, 'monitorDistInches'), settings.monitorDistInches = 36; end
if ~isfield(settings, 'monitorSizeInches'), settings.monitorSizeInches = [20.889 11.75]; end
% calculate pixels per degree.
settings.monitorPixelsPerInch = settings.screenSize(3) / settings.monitorDistInches(1);
settings.monitorPixelsPerDegree = tan(deg2rad(1)) * settings.monitorDistInches ...
    * settings.monitorPixelsPerInch;
% Interface settings - which keys do what
if ~isfield(settings, 'keyGo'), settings.keyGo = 'space'; end
if ~isfield(settings, 'keyGoName'), settings.keyGoName = 'the space bar'; end
if ~isfield(settings, 'keyLeft'), settings.keyLeft = 'LeftArrow'; end
if ~isfield(settings, 'keyLeftName'), settings.keyLeftName = 'the left arrow key'; end
if ~isfield(settings, 'keyRight'), settings.keyRight = 'RightArrow'; end
if ~isfield(settings, 'keyRightName'), settings.keyRightName = 'the right arrow key'; end
if ~isfield(settings, 'keyExit'), settings.keyExit = 'escape'; end
if ~isfield(settings, 'keyExitName'), settings.keyExitName = 'ESCAPE'; end

end