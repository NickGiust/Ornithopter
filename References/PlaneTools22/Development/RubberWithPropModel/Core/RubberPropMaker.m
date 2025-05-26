clear;
clc;

%% Import data
propNames = [];
source_dir = 'D:\Extracurricular\ADT\PlaneTools22\Rubber\Development\UIUC-propDB\volume-1\data';
propName = 'apce_11x10';
source_files = dir(fullfile(source_dir, [propName,'_*']));
opts = detectImportOptions(fullfile(source_dir, source_files(1).name));
P = [];
propNo = 0;
%propMaker = ["apc29ff", "apce", "apcsf", "apcsp", "grcp", "grsn", "gwsdd", "gwssf",...
%    "kyosho", "ma", "mae", "magf", "mas"];

propMaker = ["apce"];

for maker = 1:length(propMaker)
    for diameter = 6:1:20
        for pitch = 3:1:12
            try
                propName = [char(propMaker(maker)) '_' num2str(diameter) 'x' num2str(pitch) ];
                source_files = dir(fullfile(source_dir, [propName,'_*']));
                data = [];

                for i = 1:length(source_files)
                    if ~contains(source_files(i).name, 'geom') && ~contains(source_files(i).name, 'static')
                        data = [data; readmatrix(fullfile(source_dir, source_files(i).name), opts)];
                    end
                end

                if length(data) > 0
                    propNo = propNo + 1;
                    P = [P, atan(pitch/diameter)];
                end

                if propNo == 1
                    CP = [data(:,1), data(:,3)];
                    CT = [data(:,1), data(:,2)];
                elseif propNo ~= 0
                    CP = [CP, zeros(length(CP),1); data(:,1), zeros(length(data),propNo-1), data(:,3)];
                    CT = [CT, zeros(length(CT),1); data(:,1), zeros(length(data),propNo-1), data(:,2)];
                end
            end
        end
    end
end

CP = UniquifyDimension(CP,1);
CT = UniquifyDimension(CT,1);

CP((CP(:,1) < 0.3),:) = [];
CT((CT(:,1) < 0.3),:) = [];

J = CP(:,1);
CP = [P; CP(:,2:end)];
CT = [P; CT(:,2:end)];

CP = UniquifyDimension(CP,2);
CT = UniquifyDimension(CT,2);

P = CP(1,:);
CP = CP(2:end,:);
CT = CT(2:end,:);

CP((CP == 0)) = NaN;
CT((CT == 0)) = NaN;

[xData, yData, zData] = prepareSurfaceData( P, J, CP );
% Set up fittype and options.
ft = fittype( 'poly32' );
% Fit model to data.
[CPfit, CPgof] = fit( [xData, yData], zData, ft );
CPv = coeffvalues(CPfit)
% Plot fit with data.
figure(1);
h = plot( CPfit, [xData, yData], zData );
legend( h, 'C_{P,fit}', 'C_{P,experimental}', 'Location', 'NorthEast' );
% Label axes
xlabel('Pitch angle, \gamma [rad]')
ylabel('Advance ratio, J')
zlabel('Propeller Coefficient of Power, C_P')
grid on

[xData, yData, zData] = prepareSurfaceData( P, J, CT );
% Set up fittype and options.
ft = fittype( 'poly32' );
% Fit model to data.
[CTfit, CTgof] = fit( [xData, yData], zData, ft );
CTv = coeffvalues(CTfit)
% Plot fit with data.
figure(2);
h = plot( CTfit, [xData, yData], zData );
legend( h, 'C_{P,fit}', 'C_{P,experimental}', 'Location', 'NorthEast' );
% Label axes
xlabel('Pitch angle, \gamma [rad]')
ylabel('Advance ratio, J')
zlabel('Propeller Coefficient of Thrust, C_T')
grid on