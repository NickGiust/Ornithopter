clc;
clear;

%% Import data
propName = 'apc_20x12';
diameter = 20;

CP = [];
CT = [];

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
load('propeller.mat');
save('propeller.mat');