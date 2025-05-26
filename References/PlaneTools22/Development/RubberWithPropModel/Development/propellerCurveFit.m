clc;
clear;

%% Import data
diameter = 17;
pitch = 12;
propeller = ['apce_' int2str(diameter) 'x' int2str(pitch) '_rd'];
source_dir = 'UIUC-propDB\volume-1\data';
source_files = dir(fullfile(source_dir, [propeller,'*']));
data = [];
for i = 1:length(source_files)
    opts = detectImportOptions(fullfile(source_dir, source_files(i).name));
    data = [data; readmatrix(fullfile(source_dir, source_files(i).name), opts)];
end
CP = [data(:,1), data(:,3)];
CT = [data(:,1), data(:,2)];

%% Power Coefficient
CP = sortrows(CP);

% Values with J > .25 seem good
CP((CP(:,1) < 0.25),:) = [];

% Plot them
scatter(CP(:,1), CP(:,2));
hold on;

% Remove NaN rows from data
CP(any(isnan(CP), 2), :) = [];

% Fit both datasets
[CPfit, gof] = fit(CP(:,1),CP(:,2),'poly2');
plot(CPfit);
CPrsquare = gof.rsquare;
EFG = coeffvalues(CPfit);

title(['Aeronaut ' int2str(diameter) 'x' int2str(pitch) ' - Power']);
xlabel('Advance Ratio, \itJ');
ylabel('Power Coefficient, \itC_P');
hold off;

%% Thrust Coefficient
CT = sortrows(CT);

% Values with J > .25 seem good
CT((CT(:,1) < 0.25),:) = [];

% Plot them
figure(2);
scatter(CT(:,1), CT(:,2));
hold on;

% Remove NaN rows from data
CT(any(isnan(CT), 2), :) = [];

% Fit both datasets
[CTfit, gof] = fit(CT(:,1),CT(:,2),'poly2');
plot(CTfit);
CTrsquare = gof.rsquare;
ABC = coeffvalues(CTfit);

title([int2str(diameter) 'x' int2str(pitch) ' - Thrust']);
xlabel('Advance Ratio, \itJ');
ylabel('Thrust Coefficient, \itC_T');
hold off;

propvals = [ABC EFG]