clc;
clear;

%% Import data
propNames = [];
source_dir = 'D:\Extracurricular\ADT\PlaneTools22\Development\UIUC-propDB\volume-1\data';
propName = 'apce_11x10';
source_files = dir(fullfile(source_dir, [propName,'_*']));
opts = detectImportOptions(fullfile(source_dir, source_files(1).name));

for diameter = 6:1:20
    for pitch = 3:1:12
        try
            propName = ['zin_' num2str(diameter) 'x' num2str(pitch) ];
            source_files = dir(fullfile(source_dir, [propName,'_*']));
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
            
            propNames = [propNames ' ' propName];
        end
    end
end

eval(['clearvars -except ' propNames]);
load('propeller.mat');
save('propeller.mat');