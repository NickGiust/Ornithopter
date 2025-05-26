function [plane] = definePlane()
%DEFINEPLANE Summary of this function goes here
%   Detailed explanation goes here
plane = [];

%% WING
%=======================================================================
plane.b = 1.5; % Span [m]
plane.AR = 5.5;
plane.c = plane.b / plane.AR; % Avg. Chord [m]
plane.S = plane.b * plane.c;
plane.wingControlFraction = 0.3;
plane.tcRatio = 0.2;
plane.taperRatio = 1;
plane.croot = 2*plane.c / (1 + plane.taperRatio);
plane.ctip = plane.croot*plane.taperRatio;
% plane.airfoilType = 'e560';
plane.airfoilCoordinatesFile = 'E1098t12c2.csv';
plane.CLmax = 1.9;
plane.wingThickness = plane.c*plane.tcRatio;
plane.wingThicknessLocation = 0.24;
plane.CLgroundRoll = 0.5;
plane.hWinglet = 0;
plane.mWinglets = 0;
plane.hWing = 0.31;
plane.nStruct = 10;
plane.webMaterial = "plywood";
plane.wingBuildMethod = "builtup";
% plane.wingBuildMethod = "foamcore";
plane.controlSurfMat = "balsa";
% plane.controlSurfMat = "foam";
plane.wingMat = "carbon"; % from material.mat

%% FUSELAGE
%=======================================================================
plane.lFuse = 1; % actually calculated in TradeStudy.m
plane.wFuse = 0.115; % width of 1 patient (1" + 11/16") + 1 emt (1.5") (plus some fudge room)
plane.hFuse = 0.13; % minimum 3.5 inch (height of crew). Chose 5 inches for safety
plane.fuseMat = 'carbon';
plane.fuseBuildMethod = "monocoque"; % builtup/monocoque

%% TAIL
%=======================================================================
plane.Vh = 0.7;
plane.Vv = 0.075;
% bh = 0.4*b;
plane.ARh = 3.5;
plane.lh = 0.6*plane.b;
% Vh = Sh*lh / (S*c)
plane.Sh = plane.Vh*plane.S*plane.c/plane.lh;
plane.bh = sqrt(plane.ARh*plane.Sh);
plane.ch = plane.Sh./plane.bh;
plane.lv = plane.lh;
% Vv = Sv*lv / (S*b)
plane.Sv = plane.Vv*plane.S*plane.b/plane.lv;
plane.cv = plane.ch;
plane.bv = plane.Sv/plane.cv;
plane.ARv = plane.bv.^2./plane.Sv;
% Tail
plane.cTail = plane.ch;
plane.bTail = plane.bh;
plane.hTail = plane.bv;
plane.nTails = 1;
plane.tailMat = "builtup";
% plane.bTail = 0.47;
% plane.hTail = 0.27;
% plane.cTail = 0.15;

%% PROPULSION
%=======================================================================
plane.motorType = 'Scorpion_A_4220_540';
% plane.motorType = 'Scorpion_HKIII_5025_520KV_F3S';
warning('off')
load('motor.mat', plane.motorType);
motor = eval(plane.motorType);
plane.Kv = motor.Kv;
plane.motorMaxPower = motor.maxPower;
Rmotor = motor.R;
plane.motorLength = motor.length;
plane.motorDiam = motor.diam;
plane.nMotors = 2;
plane.ESCMaxCurrent = 260;
plane.ESCContinuousCurrent = 110;
Resc = .0005;
Rwire = .01;

%% LANDING GEAR
%=======================================================================
plane.lgType = 'bow';
plane.lgMat = 'carbon';
plane.rr = 0.03; % rolling resistance

%% MISSION 2
%=======================================================================
% Propeller
diameter = 11;
pitch = 7;
% diameter = 14;
% pitch = 8;
propellerMAT = matfile('myFile2.mat');
allProps = propellerMAT.allProps;
diamIndices = find(allProps(:,1)==diameter);
pitchIndices = find(allProps(:,2)==pitch);
[index,~]=intersect(diamIndices,pitchIndices);

plane.D2 = diameter*0.0254;
plane.P2 = pitch;
plane.A2 = allProps(index,3);
plane.B2 = allProps(index,4);
plane.C2 = allProps(index,5);
plane.E2 = allProps(index,6);
plane.F2 = allProps(index,7);
plane.G2 = allProps(index,8);

% Battery
plane.batteryType2 = 'Thunder_3300_8s';
load('battery.mat', plane.batteryType2);
battery2 = eval(plane.batteryType2);

plane.bat2maxCurrent = battery2.Imax;
plane.bat2R = battery2.R;

% Total energy
plane.nSeries2 = battery2.nSeries;
plane.nParallel2 = 1;
plane.bat2capacity = battery2.capacity*plane.nParallel2;
plane.Eb2 = plane.bat2capacity*plane.nSeries2*3.7*3600;

% Total resistance
plane.Rt2 = Resc + Rwire + Rmotor + (plane.bat2R/plane.nParallel2);

% Payload
plane.nPayloads2 = 1;
plane.massPayloads2 = 0.35; % overwritten
plane.DApayload2 = 0;
mPerson = 0.04; % [lbs -> kg] (mass of 25 dolls) / 25
mPatient = 1.52 / 5 / 2.205; % [lbs -> kg] (mass of 5 dolls) / 5
plane.mPayloadExtras2 = (4 * mPerson) + mPatient + 0.1; % [kg], mass of 2 crew, 2 emt, 1 patient, 1 gourney

plane.seed_m2 = 6;

% Flight Performance
plane.Vmax2 = 80 / 2.237; % [mph -> m/s], max speed plane should be flying at

%% MISSION 3
%=======================================================================
% Propeller
diameter = 11;
pitch = 6;
% diameter = 14;
% pitch = 8;
propellerMAT = matfile('myFile2.mat');
allProps = propellerMAT.allProps;
diamIndices = find(allProps(:,1)==diameter);
pitchIndices = find(allProps(:,2)==pitch);
[index,~]=intersect(diamIndices,pitchIndices);

plane.D3 = diameter*0.0254;
plane.P3 = pitch;
plane.A3 = allProps(index,3);
plane.B3 = allProps(index,4);
plane.C3 = allProps(index,5);
plane.E3 = allProps(index,6);
plane.F3 = allProps(index,7);
plane.G3 = allProps(index,8);

% Battery
plane.batteryType3 = 'Thunder_3300_8s'; 
load('battery.mat', plane.batteryType3);
battery3 = eval(plane.batteryType3);

plane.bat3maxCurrent = battery3.Imax;
plane.bat3R = battery3.R;

% Total energy
plane.nSeries3 = battery3.nSeries;
plane.nParallel3 = 1;
plane.bat3capacity = battery3.capacity*plane.nParallel3;
plane.Eb3 = plane.bat3capacity*plane.nSeries3*3.7*3600;

% Total resistance
plane.Rt3 = Resc + Rwire + Rmotor + (plane.bat3R/plane.nParallel3);

% Payload
plane.nPayloads3 = 3; % overwritten
plane.massPayloads3 = mPerson; % [kg]
plane.DApayload3 = 0;
plane.mPayloadExtras3 = (2 * mPerson); % [kg], mass of 2 crew

% Flight Performance
plane.Vmax3 = 80 / 2.237; % [mph -> m/s], max speed plane should be flying at

%% GROUND MISSION
%=======================================================================


end

