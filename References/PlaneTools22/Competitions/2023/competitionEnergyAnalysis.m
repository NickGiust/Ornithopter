
clear; close all; clc;


%% limits as given by the rules

% bLimit = 4*0.3048; % 4 ft span limit
% b = bLimit;
energyLimitWh = 100; % 100 Wh energy limit
energyLimitJ = energyLimitWh*3600; % energy limit in Joules
TOFLLimit = 60*0.3048; % 60 ft takeoff field length limit
distanceHalfStraight = 500*0.3048; % full straightaways are 1000 ft long
timeLimit = 600;
nLapsM3 = 3;


%% constants and efficiencies

g = 9.81;
Crr = 0.05; % coefficient of rolling resistance
% find values here: https://www.engineeringtoolbox.com/rolling-friction-resistance-d_1303.html
rho = 1.1303; % tucson at ~800 m altitude, using stdatm function from AME 261
dens = rho;
turnReactionTime = 0.5;
% e = 0.9; % formula used to find e

etaESC = 0.95;
etaBattery = 0.97;
etaMotor = 0.94;
etaPropeller = 0.8;
etaTotal = etaESC*etaBattery*etaMotor*etaPropeller;

batteryEndSoC = 0.1; % discharge battery down to 10% SoC


%% parameter sweep list



bList = 0.5:0.25:2.5; % span will always be the maximum to reduce induced drag
ARList = 4.5:0.25:7; % aspect ratio between 3 and 15, seems reasonable
CLmaxList = 2.0; % maximum lift coefficient between 1.4 and 2.2
CLCruiseFullList = 0.1:0.05:0.8; % average cruise lift coefficient while heavy
CLCruiseEmptyList = 0.1:0.05:0.4; % average cruise lift coefficient after dumping water
TStaticList = 70; % static thrust, N ( = 1 kg to 20 kg)
% TDynamicShareFullList = 0.05:0.05:0.2; % drag in cruise while heavy
% TDynamicShareEmptyList = 0.05:0.05:0.2; % drag in cruise after dumping water
nTurnFullList = 4:0.5:10; % number of G's pulled in turns while full, assuming level
nTurnEmptyList = 3:0.5:9; % number of G's pulled in turns while empty, assuming level
VWaterList = 0.5:0.25:3; % volume of water payload, Liters
numLapsList = 12:2:32; % total number of laps flown
buildTechnologyFactorList = 0.9; % multiplier on aircraft empty weight




% Since T = D = 0.5*rho*U^2*S*(CL^2/(pi*e*AR) + CDO)
% density is fixed, speed is unknown, S is swept over, CL is unknown
% e estimated as a constant value, AR is swept over, and CDO = f(geometry)
% which is known (assuming speed is estimated for Re calculation)
% In summary, thrust/drag, speed, and lift coefficient are unknown.
% Therefore, two of thrust/speed/lift must be set to have a determinate
% (steady state) cruise condition
% We will choose to set the cruise lift coefficient


% create a unique plane for every permutation
[AR, CLmax, CLCruiseFull, CLCruiseEmpty, TStatic, nTurnFull, nTurnEmpty,...
    VWater, numLaps, buildTechnologyFactor] = ...
    ndgrid(ARList, CLmaxList, CLCruiseFullList, CLCruiseEmptyList,...
    TStaticList, nTurnFullList, nTurnEmptyList, VWaterList, numLapsList,...
    buildTechnologyFactorList);

%[AR, CLmax, CLCruiseFull, CLCruiseEmpty, TStatic, TDynamicShareFull,...
%     TDynamicShareEmpty, nTurnFull, nTurnEmpty, VWater, numLapsList] = ...
%     ndgrid(ARList, CLmaxList, CLCruiseFullList, CLCruiseEmptyList,...
%     TStaticList, TDynamicShareFullList, TDynamicShareEmptyList,...
%     nTurnFullList, nTurnEmptyList, VWaterList, numLapsList);

numLapsFull = ceil(numLaps/2); % at least half of laps must be with water
numLapsEmpty = numLaps - numLapsFull; % at max half of laps are empty

numPlanes = prod(size(AR), 'all');

fprintf("Running %.1f million planes\n\n", numPlanes/(1e6));

S = b^2./AR;
c = b./AR;


%% calculate fuselage dimensions

VFuseLiters = 1.5.*VWater + 3; % estimating 2 L for battery, avionics, payload servos
    % also with 50% extra volume to account for the tailcone

VFuse = VFuseLiters/1000; % convert fuselage volume in liters to m^3

finenessFuse = 6; % arbitrarily choosing fuselage fineness ratio

% Volume of ellipsoid is pi*A*B*C/6, where A, B, and C are the principal
% diameters (see https://en.wikipedia.org/wiki/Ellipsoid)
% Here we let B = C since we approximate the fuselage to have the same
% height and width. Therefore VFuse = pi*lFuse*dFuse^2/6 where lFuse is the
% fuselage length and dFuse is the widest fuselage diameter.
% Since finenessFuse = lFuse/dFuse, lFuse = dFuse*finenessFuse.
% Therefore, VFuse = pi*finenessFuse*dFuse^3/6. Since volume and fineness
% ratio are already determined, solve for dFuse such that:
% dFuse = (6*VFuse / (pi*finenessFuse))^(1/3)

dFuse = (6*VFuse ./ (pi*finenessFuse)) .^ (1/3);
lFuse = dFuse*finenessFuse;











%% DRAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taken from PlaneTools22 dragBuildup function



t = c.*0.16; % assuming e560 with 16.1% t/c @ 24% chord
tLocation = 0.24;
% http://airfoiltools.com/airfoil/details?airfoil=e560-il

tTail = cTail*.1;


% Air
Vapprox = 25; % using this speed in m/s to find Re for CD0 calculations
tempC = 18;
tempK = tempC + 273.15;
a = sqrt(1.4*287*tempK);
Mach = Vapprox/a;



%% Reynolds Numbers
% Wing
ReWing = (dens*Vapprox.*c)./1.8e-5;

% Fuselage
ReFuse = (dens*Vapprox.*lFuse)./1.8e-5;

% Tail 
ReTail = (dens*Vapprox*cTail)/1.8e-5;

%% Oswald efficiency
% Initial rectangular planform
f = .0524*(1^4) - 0.15*(1^3) + .1659*(1^2) - .0706*1 + .0119; % assuming taper ratio of 1 here
e = 1./(1 + (f.*AR));


%% Wing
% calculate laminar and turbulent skin friction coefficients
% Prandtl
CfWingSmoothLaminar = 1.328./sqrt(ReWing);
% Roskam
CfWingSmoothTurbulent = 0.455 ./ (log10(ReWing).^2.58 .* (1+0.144.*Mach^2).^0.65);

proportionLaminar = 0.6;
CfWingSmooth = proportionLaminar*CfWingSmoothLaminar +...
    (1-proportionLaminar)*CfWingSmoothTurbulent;

% L' = 1.2 for (t/c)max at x_t >= 0.3c
% L' = 2.0 for (t/c)max at x_t < 0.3c
% ex. BA527ls has (t/c)max at x_t = 34%c, so L' = 1.2
% ex. e560 has (t/c)max at x_t = 24%c, so L' = 2.0

% the location of max thickness is aft of 30% chord
if tLocation >= 0.3
    Lprime = 1.2;
else
    % max thickness is forward of 30% chord
    Lprime = 2.0;
end

R_wf = 1.1; % Roskam page 24, taking Re_fus ~ 2-3e6 and M<0.25
R_ls = 1.07; % Roskam page 24, based on zero sweep and M<0.25

SwetWing = 2*S*1.05;

FFwing = (1 + Lprime*(t./c) + 100*(t./c).^4);

DAwing = R_wf*R_ls .* CfWingSmooth .* FFwing .* SwetWing;




%% Tail
% calculate laminar and turbulent skin friction coefficients
% Prandtl
CfTailSmoothLaminar = 1.328/sqrt(ReTail);
% Roskam
CfTailSmoothTurbulent = 0.455 ./ (log10(ReTail).^2.58 .* (1+0.144*Mach^2).^0.65);


CfTailSmooth = proportionLaminar*CfTailSmoothLaminar +...
    (1-proportionLaminar)*CfTailSmoothTurbulent;


% Form factor, Roskam (1.6 is L' for NACA 0010 with (t/c)max at 30%c
FFtail = 1 + (1.6.*(tTail./cTail)) + (100.*((tTail./cTail).^4));

% Total tail drag area;
DAtail = 2.05.*(cTail.*bTail.*R_ls + cTail.*hTail).*CfTailSmooth.*FFtail;
% accounting for lifting surface correction, Roskam page 24


%% Fuselage

% Roskam
CfFuseSmoothTurbulent = 0.455 ./ (log10(ReFuse).^2.58 * (1+0.144*Mach.^2).^0.65);

CfFuse = CfFuseSmoothTurbulent; % assuming fuselage has entirely turbulent boundary layer

% form factor, Roskam page 44
FFfuse = (1 + 60./(lFuse./dFuse).^3 + 0.0025.*(lFuse./dFuse));

% drag area, Roskam page 44, ignoring C_D_b_fus (small contribution)
DAfuse = R_wf.*CfFuse .* FFfuse .* SwetFuse;


%% Landing Gear

CDlg = .1; % applicable for bow
Alg = hWing*(hWing/10);
DAlg = CDlg*Alg;
DAlg = (AR.*0 + 1)*DAlg;

%% Total
DA = DAwing + DAfuse + DAtail + DAlg;

CD0 = DA./S;








