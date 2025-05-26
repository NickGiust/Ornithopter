function plane = ReadPlane(planeName, directory)
%% All missions
% Text file import as table
addpath([directory '\Library']);
opts = detectImportOptions(fullfile([directory '\Planes\' planeName '.txt']));
opts.VariableOptions(1,2) = opts.VariableOptions(1,1);
txtParams = readtable(fullfile([directory '\Planes\' planeName '.txt']), opts);
params = str2double(txtParams{:,2});

% Name
plane = [];
plane.name = planeName;

% Wing
plane.b = params(1);
plane.cTip = params(2);
plane.cRoot = params(3);
plane.c = (plane.cTip + plane.cRoot)/2;
plane.t = params(4);

plane.airfoilType = char(txtParams{5,2});
load('airfoil.mat', plane.airfoilType);
airfoil = eval(plane.airfoilType);
plane.CLmax = airfoil.CLmax2D;
plane.wingThickness = airfoil.thicknessRatio*plane.c;
plane.CLgroundRoll = params(6);

plane.DAwinglets = params(7);
plane.mWinglets = params(8);

plane.hWing = params(9);
plane.nStruct = params(10);
plane.wingMat = char(txtParams{11,2});

% Fuselage
plane.lFuse = params(14);
plane.wFuse = params(15);
plane.hFuse = params(16);
plane.fuseMat = char(txtParams{17,2});

% Tail
plane.bTail = params(20);
plane.hTail = params(21);
plane.cTail = params(22);
plane.nTails = params(23);
plane.tailMat = char(txtParams{24,2});

% Motor
plane.motorType = char(txtParams{27,2});
load('motor.mat', plane.motorType);
motor = eval(plane.motorType);

plane.Kv = motor.Kv;
plane.motorMaxPower = motor.maxPower;
Rmotor = motor.R;
plane.nMotors = params(28);

% ESC
plane.ESCMaxCurrent = params(29);
Resc = params(30);

% Wire
Rwire = .01;

% Landing gear
plane.lgType = char(txtParams{34,2});
plane.lgMat = char(txtParams{35,2});
plane.rr = params(36);

%% Mission 2
% Propeller
plane.propType2 = char(txtParams{40,2});
load('propeller.mat', plane.propType2);
propeller = eval(plane.propType2);

plane.D2 = propeller.D;
plane.A2 = propeller.A;
plane.B2 = propeller.B;
plane.C2 = propeller.C;
plane.E2 = propeller.E;
plane.F2 = propeller.F;
plane.G2 = propeller.G;

% Battery
plane.batteryType2 = char(txtParams{42,2});
load('battery.mat', plane.batteryType2);
battery2 = eval(plane.batteryType2);

capacity2 = battery2.capacity;
plane.bat2maxCurrent = battery2.Imax;
Rbat2 = battery2.R;

% Total energy
plane.nSeries2 = params(43);
plane.nParallel2 = params(44);
plane.Eb2 = capacity2*3.7*plane.nSeries2*plane.nParallel2*3600;

% Total resistance
plane.Rt2 = Resc + Rwire + Rmotor + Rbat2;

% Payload
plane.nPayloads2 = params(48);
plane.mPayload2 = params(46);
plane.DApayload2 = params(47);

%% Mission 3
% Propeller
plane.propType3 = char(txtParams{51,2});
load('propeller.mat', plane.propType3);
propeller = eval(plane.propType3);

plane.D3 = propeller.D;
plane.A3 = propeller.A;
plane.B3 = propeller.B;
plane.C3 = propeller.C;
plane.E3 = propeller.E;
plane.F3 = propeller.F;
plane.G3 = propeller.G;

% Battery
plane.batteryType3 = char(txtParams{53,2});
load('battery.mat', plane.batteryType3);
battery3 = eval(plane.batteryType3);

capacity3 = battery3.capacity;
plane.bat3maxCurrent = battery3.Imax;
Rbat3 = battery3.R;

% Total energy
plane.nSeries3 = params(54);
plane.nParallel3 = params(55);
plane.Eb3 = capacity3*3.7*plane.nSeries3*plane.nParallel3*3600;

% Total resistance
plane.Rt3 = Resc + Rwire + Rmotor + Rbat3;

% Payload
plane.nPayloads3 = params(59);
plane.mPayload3 = params(57);
plane.DApayload3 = params(58);

end