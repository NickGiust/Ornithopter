source_dir = 'D:\ADT\PlaneTools22\Development\UIUC-propDB\volume-3\data';
propName = 'ancf_15x10';
diameter = 15;
source_files = dir(fullfile(source_dir, [propName,'_*']));!excel.exe &
opts = detectImportOptions(fullfile(source_dir, source_files(1).name));

data = [];

for i = 1:length(source_files)
    if ~contains(source_files(i).name, 'geom') && ~contains(source_files(i).name, 'static')
        data = [data; readmatrix(fullfile(source_dir, source_files(i).name), opts)];
    end
end

CP = [data(:,1), data(:,3)];
CT = [data(:,1), data(:,2)];

%% Power Coefficient
CP = sortrows(CP);

% Values with J > .25 seem good
CP((CP(:,1) < 0.25),:) = [];

% Remove NaN rows from data
CP(any(isnan(CP), 2), :) = [];

% Fit both datasets
CPfit = fit(CP(:,1),CP(:,2),'poly2');

EFG = coeffvalues(CPfit);

%% Thrust Coefficient
CT = sortrows(CT);

% Values with J > .25 seem good
CT((CT(:,1) < 0.25),:) = [];

% Remove NaN rows from data
CT(any(isnan(CT), 2), :) = [];

% Fit both datasets
CTfit = fit(CT(:,1),CT(:,2),'poly2');
ABC = coeffvalues(CTfit);

propeller.D = diameter*.0254;
propeller.A = ABC(1);
propeller.B = ABC(2);
propeller.C = ABC(3);
propeller.E = EFG(1);
propeller.F = EFG(2);
propeller.G = EFG(3);
if contains(propName, '.')
    propName = replace(propName, '.', 'p');
end

eval([propName ' =  propeller; clear propeller'])

eval(['clearvars -except ' propName]);