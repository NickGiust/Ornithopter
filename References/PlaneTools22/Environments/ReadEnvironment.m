function environment = ReadEnvironment(envName, directory)
% global slash
% slash
currentDirectory = split(pwd, 'PlaneTools22');
currentDirectoryCharacterVector = char(currentDirectory(1,:));
backslashLocations = strfind(currentDirectoryCharacterVector, '\');
if isempty(backslashLocations)
    slash = '/'; % forward slash for Apple Mac directories
else
    slash = '\'; % back slash for Windows directories
end
% Text file import as table
opts = detectImportOptions(fullfile([directory slash 'Environments' slash envName '.txt']));
opts.VariableOptions(1,2) = opts.VariableOptions(1,1);
txtParams = readtable(fullfile([directory slash 'Environments' slash envName '.txt']), opts);
params = str2double(txtParams{:,2});

%% Environment parameters
environment = [];
windDirectionDegrees = params(3);
runwayDirectionUpwindDegrees = params(4);
environment.windSpeed = params(2);
environment.dens = params(1);

environment.windDirection = deg2rad(windDirectionDegrees);
environment.runwayDirectionUpwind = deg2rad(runwayDirectionUpwindDegrees);