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
if ~isfield(settings, 'ptbSkipSyncTests'), settings.ptbSkipSyncTests = false; end
% Interface settings - which keys do what
if ~isfield(settings, 'keyGo'), settings.keyGo = 'space'; end
if ~isfield(settings, 'keyGoName'), settings.keyGoName = 'the space bar'; end
if ~isfield(settings, 'keyLeft'), settings.keyLeft = 'left'; end
if ~isfield(settings, 'keyLeftName'), settings.keyLeftName = 'LEFT'; end
if ~isfield(settings, 'keyRight'), settings.keyRight = 'right'; end
if ~isfield(settings, 'keyRightName'), settings.keyRightName = 'RIGHT'; end
if ~isfield(settings, 'keyExit'), settings.keyExit = 'escape'; end
if ~isfield(settings, 'keyExitName'), settings.keyExitName = 'ESCAPE'; end

end