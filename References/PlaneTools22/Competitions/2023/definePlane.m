function [plane] = definePlane()
%DEFINEPLANE Summary of this function goes here
%   Detailed explanation goes here
plane = [];

%% WING
%=======================================================================
plane.b = 1.25; % Span [m]
% plane.cTip = 0.21; % Tip chord [m]
% plane.cRoot = 0.21; % Root chord [m]
plane.AR = 5.5;
plane.c = plane.b / plane.AR; % Avg. Chord [m]
plane.S = plane.b * plane.c;
plane.wingControlFraction = 0.5;
plane.tcRatio = 0.2;
plane.taperRatio = 0.8;
plane.croot = 2*plane.c / (1 + plane.taperRatio);
plane.ctip = plane.croot*plane.taperRatio;
plane.airfoilType = 'e560';
plane.airfoilCoordinatesFile = 'goe226-il.csv';
plane.CLmax = 1.7;
plane.wingThickness = plane.c*plane.tcRatio;
plane.wingThicknessLocation = 0.24;
plane.CLgroundRoll = 1.2;
plane.hWinglet = 0;
plane.mWinglets = 0;
plane.hWing = 0.31;
plane.nStruct = 10;
% plane.wingMat = 'foamfiberglass';
plane.webMaterial = "plywood";
plane.wingBuildMethod = "builtup"; % "builtup" OR "foamcore"
plane.controlSurfMat = "balsa";
plane.numberOfTubes = 2;



%% FUSELAGE
%=======================================================================
plane.lFuse = 0.65;
plane.wFuse = 0.15;
plane.hFuse = 0.15;
plane.fuseMat = 'carbon';

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
plane.motorType = 'Neu_Motors_1915_1Y';
load('motor.mat', plane.motorType);
motor = eval(plane.motorType);
plane.Kv = motor.Kv;
plane.motorMaxPower = motor.maxPower;
Rmotor = motor.R;
plane.nMotors = 1;
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
diameter = 13;
pitch = 10;
propellerMAT = matfile('myFile2.mat');
allProps = propellerMAT.allProps;
diamIndices = find(allProps(:,1)==diameter);
pitchIndices = find(allProps(:,2)==pitch);
[index,pos]=intersect(diamIndices,pitchIndices);

plane.D2 = diameter*0.0254;
plane.P2 = pitch;
plane.A2 = allProps(index,3);
plane.B2 = allProps(index,4);
plane.C2 = allProps(index,5);
plane.E2 = allProps(index,6);
plane.F2 = allProps(index,7);
plane.G2 = allProps(index,8);

% Battery
plane.batteryType2 = 'Thunder_LiPo_2700';
load('battery.mat', plane.batteryType2);
battery2 = eval(plane.batteryType2);

plane.bat2maxCurrent = battery2.Imax;
plane.bat2R = battery2.R;

% Total energy
plane.nSeries2 = 10;
plane.nParallel2 = 1;
plane.bat2capacity = battery2.capacity*plane.nParallel2;
plane.Eb2 = plane.bat2capacity*plane.nSeries2*3.7*3600;

% Total resistance
plane.Rt2 = Resc + Rwire + Rmotor + (plane.bat2R/plane.nParallel2);

% Payload
plane.nPayloads2 = 1;
plane.payload2Fraction = 0.35;
plane.DApayload2 = 0;

plane.seed_m2 = 6;

%% MISSION 3
%=======================================================================
% Propeller
diameter = 14;
pitch = 12;
propellerMAT = matfile('myFile2.mat');
allProps = propellerMAT.allProps;
diamIndices = find(allProps(:,1)==diameter);
pitchIndices = find(allProps(:,2)==pitch);
[index,pos]=intersect(diamIndices,pitchIndices);

plane.D3 = diameter*0.0254;
plane.P3 = pitch;
plane.A3 = allProps(index,3);
plane.B3 = allProps(index,4);
plane.C3 = allProps(index,5);
plane.E3 = allProps(index,6);
plane.F3 = allProps(index,7);
plane.G3 = allProps(index,8);

% Battery
plane.batteryType3 = 'Thunder_LiPo_2700'; 
load('battery.mat', plane.batteryType3);
battery3 = eval(plane.batteryType3);

plane.bat3maxCurrent = battery3.Imax;
plane.bat3R = battery3.R;

% Total energy
plane.nSeries3 = 10;
plane.nParallel3 = 1;
plane.bat3capacity = battery3.capacity*plane.nParallel3;
plane.Eb3 = plane.bat3capacity*plane.nSeries3*3.7*3600;

% Total resistance
plane.Rt3 = Resc + Rwire + Rmotor + (plane.bat3R/plane.nParallel3);

% Payload
plane.nPayloads3 = 1;
plane.lengthPayload3 = 0.8;
plane.DApayload3 = (1.17*(plane.lengthPayload3*0.0214))/plane.S;

%% GROUND MISSION
%=======================================================================
plane.GMloadingMultiplier = 10; % multiplier on M2 mass

end

