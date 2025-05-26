function [massResults, plane] = massBuildup(plane, mission)

%% PARAMETERS

if plane.wingBuildMethod.isequal("builtup")
    [mWing, plane] = calcMass_BuiltupWing(plane);
    plane.mWing = mWing;
elseif plane.wingBuildMethod.isequal("foamcore")
    [mWing, plane] = calcMass_FoamcoreWing(plane);
    plane.mWing = mWing;
end

% Winglet
density_balsa = 160; %kg/m3
cWinglet = 0.5*(plane.ctip+0.5*plane.ctip);
Swinglet = cWinglet*plane.hWinglet;
wingletAvgThickness = 0.08 * cWinglet;
mWinglet = wingletAvgThickness*Swinglet*2*density_balsa;
plane.mWinglet = mWinglet;

% Fuselage
mFuse = calcMass_Fuselage(plane);
plane.mFuse = mFuse;

% Conventional Tail / T Tail (Built-up)

% builtupDensity = mWing / (plane.b * plane.c); % kg/m2
% mVStab = builtupDensity * plane.cTail * plane.hTail;
% mHStab = builtupDensity * plane.cTail * plane.bTail;
% mTail = mVStab + mHStab;

tailAvgThickness = 0.05 * plane.cTail;
mTail = plane.cTail*(plane.bTail + plane.hTail)*tailAvgThickness*density_balsa;
plane.mTail = mTail;
% V Tail (Built-up)

% Landing gear - Sarahhh!!!
lgType = plane.lgType;
hWing = plane.hWing;

if strcmp(lgType, 'strut')
    Alg = (hWing^3)*(.15^2);
end

if strcmp(lgType, 'bow')
    Alg = hWing*0.2*3.5;
end

load('material.mat', plane.lgMat);
lgMaterial = eval(plane.lgMat);
daLg = lgMaterial.density;

mLg = daLg*Alg;
plane.mLg = mLg;

% Motor
load('motor.mat', plane.motorType);
motor = eval(plane.motorType);
mMotor = motor.m;
nMotors = plane.nMotors;
plane.mMotor = mMotor*nMotors;

% Motor Mount
mMotorMount = calcMass_MotorMount(plane);
plane.mMotorMount = mMotorMount;

% Battery
if mission == 2
    load('battery.mat', plane.batteryType2);
    battery = eval(plane.batteryType2);
    mBat = battery.mass*plane.nParallel2;
elseif mission == 3
    load('battery.mat', plane.batteryType3);
    battery = eval(plane.batteryType3);
    mBat = battery.mass*plane.nParallel3;
end
% mBat = 0.3; % [kg], random rough weight
plane.mBat = mBat;

% ESC
mESC = calcMass_ESC(plane, plane.ESCContinuousCurrent);
plane.mESC = mESC*nMotors;

% Avionics
mAvionics = 0.3; % kg, allotted weight for avionics
plane.mAvionics = mAvionics;

% Servo
m_tailServos = 2 * calcMass_Servos(4.8, 4.8); % 1 HS, 1 VS
m_wingServos = 4 * calcMass_Servos(4.8, 6.0); % 2 aileron, 2 flap
mServos = m_wingServos + m_tailServos;
plane.mServos = mServos;

% Wiring
mWiring = calcMass_Wiring(plane);
plane.mWiring = mWiring;

% SUM BEFORE PAYLOAD
m_emptyList = [mFuse mWing mTail mBat mMotor*nMotors mMotorMount mLg mServos mWiring mESC mAvionics]; %mWinglet];
m_empty = sum(m_emptyList);
plane.m_empty = m_empty;
plane.m_emptyList = m_emptyList;

% Payload
plane.mPayload2 = plane.nPayloads2 * plane.massPayloads2 + plane.mPayloadExtras2;
plane.mPayload3 = plane.nPayloads3 * plane.massPayloads3 + plane.mPayloadExtras3;

%% Total
plane.mContributors = ["Fuselage", "Wing", "Tail", "Battery", "Motors", "Motor Mounts", "Landing Gear", "Servos", "Wiring", "ESC", "Avionics"];

if mission == 2
    massResults.m = m_empty + plane.mPayload2;
else
    massResults.m = m_empty + plane.mPayload3;
end

end