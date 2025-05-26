function massResults = massBuildup(plane, mission)
%% PARAMETERS
% Fuselage
lFuse = plane.lFuse;
wFuse = plane.wFuse;
hFuse = plane.hFuse;

load('material.mat', plane.fuseMat);
fuseMaterial = eval(plane.fuseMat);
dFuse = fuseMaterial.density;

hollowFuse = .4;
Vfuse = (lFuse*wFuse*hFuse)*(1-hollowFuse);
mFuse = dFuse*Vfuse;

% Wings
b = plane.b;
c = plane.c;
mWinglets = plane.mWinglets;
load('material.mat', plane.wingMat);
wingMaterial = eval(plane.wingMat);
dWing = wingMaterial.density;

Vwing = b*c/60;
mWings = dWing*Vwing;

% Tail
cTail = plane.cTail;
bTail = plane.bTail;
hTail = plane.hTail;

load('material.mat', plane.tailMat);
tailMaterial = eval(plane.tailMat);
dTail = tailMaterial.density;

Vtail = cTail*(bTail+hTail)/60;
mTail = dTail*Vtail;

% Motor
load('motor.mat', plane.motorType);
motor = eval(plane.motorType);
mMotor = motor.m;

% Battery
if mission == 2
    load('battery.mat', plane.batteryType2);
    battery = eval(plane.batteryType2);
    mBat = battery.m*plane.nSeries2*plane.nParallel2;
elseif mission == 3
    load('battery.mat', plane.batteryType3);
    battery = eval(plane.batteryType3);
    mBat = battery.m*plane.nSeries3*plane.nParallel3;
end

% Landing gear
lgType = plane.lgType;
hWing = plane.hWing;

if strcmp(lgType, 'strut')
    Vlg = (hWing^3)*(.08^2);
end

load('material.mat', plane.lgMat);
lgMaterial = eval(plane.lgMat);
dLg = lgMaterial.density;

mLg = dLg*Vlg;

% Payloads
if mission == 2
    mPayloads = plane.nPayloads2*plane.mPayload2;
elseif mission == 3
    mPayloads = plane.nPayloads3*plane.mPayload3;
end

%% Total
m = [mFuse mWings mTail mBat mMotor mLg mPayloads mWinglets];
massResults.m = sum(m);
massResults.mShares = m/massResults.m;
end