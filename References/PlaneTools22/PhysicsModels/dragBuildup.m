function dragResults = dragBuildup(plane, mission)

showDragBuildup = false; % if true, displays the drag buildup every time this function is called

%% PARAMETERS
% Wing
b = plane.b;
% cavg = plane.c;
AR = plane.AR;
S = plane.S;
c = plane.c;
tcRatio = plane.tcRatio;
taperRatio = plane.taperRatio;
tLocation = plane.wingThicknessLocation;
hWinglet = plane.hWinglet;

% croot = 2*cavg / (1 + taperRatio);
% ctip = croot*taperRatio;

% load('material.mat', plane.wingMat);
% wingMaterial = eval(plane.wingMat);
% roughnessWing = wingMaterial.roughness;

% Tail
cTail = plane.cTail;
bTail = plane.bTail;
hTail = plane.hTail;
tTail = 0.1*cTail;

% load('material.mat', plane.tailMat);
% tailMaterial = eval(plane.tailMat);
% roughnessTail = tailMaterial.roughness;

% Fuselage
lFuse = plane.lFuse;
wFuse = plane.wFuse;
hFuse = plane.hFuse;

% load('material.mat', plane.fuseMat);
% fuselageMaterial = eval(plane.fuseMat);
% roughnessFuse = fuselageMaterial.roughness;

% Payloads
if mission == 2
    %% LF Radar
    DApayloads = plane.DApayload2;
elseif mission == 3
    DApayloads = plane.DApayload3;
end

% Landing gear
lgType = plane.lgType;
hWing = plane.hWing;

% Air
if mission == 2
    Vapprox = 30; % [m/s]
elseif mission == 3
    Vapprox = 40; % [m/s]
end
dens = 1.2; % [kg/m^3]
visc = 1.8e-5; % dynamic viscosity [kg/m*s]
tempC = 20; % [C]
tempK = tempC + 273.15; % [C -> K]
a = sqrt(1.4*287*tempK); % ([-] * [J/kg*K] * [K])^0.5 = ([J/kg])^0.5 = [m/s]
Mach = Vapprox/a; % [-]
ReTransition = 5e5; % [-]

% Basic equations
% c = (2/3)*croot*((taperRatio^2+taperRatio+1)/(taperRatio+1)); % unswept
% AR = b/c;
% S = b*c;

% cwinglet = 0.5*(c+0.5*c);
cwinglet = 0.75*c; % This is a simplified version of the equation above. Not sure why it's not based on winglet taper ratio
Swinglet = cwinglet*hWinglet;

%% Reynolds Numbers
% Wing
ReWing = (dens*Vapprox*c)/visc;

% Winglets
ReWinglet = (dens*Vapprox*cwinglet)/visc;

% Fuselage
ReFuse = (dens*Vapprox*lFuse)/visc;

% Tail 
ReTail = (dens*Vapprox*cTail)/visc;

%% Oswald efficiency
% Initial rectangular planform
f = .0524*(taperRatio^4) - 0.15*(taperRatio^3) + .1659*(taperRatio^2) - .0706*taperRatio + .0119;
eInitial = 1/(1 + (f*AR));

% Sweep
sweepEffect = .92; % what is this?

% Final
e = eInitial*sweepEffect;
ARmultiplier = (1 + 1.9*hWinglet/b); % Hoerner 3.9
dragResults.e = ARmultiplier*e;

%% Wing
% calculate laminar and turbulent skin friction coefficients
% Prandtl
CfWingSmoothLaminar = 1.328/sqrt(ReWing);
% Roskam
CfWingSmoothTurbulent = 0.455 / (log10(ReWing)^2.58 * (1+0.144*Mach^2)^0.65);

% flow transitions from laminar to turbulent at 500,000
if ReWing < ReTransition
    % if Reynolds number is less than 500k, fully laminar flow
    CfWingSmooth = CfWingSmoothLaminar;
else
    % Reynolds number above 500k, first part of wing in laminar flow while
    % the rest is in turbulent conditions
    CfWingSmoothProportionLaminar = ReTransition/ReWing;
    CfWingSmooth = CfWingSmoothLaminar*CfWingSmoothProportionLaminar + (1-CfWingSmoothProportionLaminar)*CfWingSmoothTurbulent;
end

% Or, with a rough surface, the following approximation from Schlichting
% https://people.csail.mit.edu/jaffer/convect/rough.pdf
% roughness is equivalent to the grit of sandpaper
% for example, roughness of 1*10^-3 = 0.001 corresponds to 1000 grit
% sandpaper
% CfWingRough = (1.89 + (1.62 * log10(c/roughnessWing)))^-2.5;
CfWingRough = 0; % NOT USED, JTS - why not???

% The greater of the two
CfWing = max([CfWingSmooth CfWingRough]);

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
% SwetWing = plane.airfoilPerimeter * b; % roughly half of the above value for some reason

FFwing = (1 + Lprime*(tcRatio) + 100*(tcRatio)^4);

DAwing = R_wf*R_ls * CfWing * FFwing * SwetWing;


%% Winglets
% calculate laminar and turbulent skin friction coefficients
% Prandtl
CfWingletSmoothLaminar = 1.328/sqrt(ReWinglet);
% Roskam
CfWingletSmoothTurbulent = 0.455 / (log10(ReWinglet)^2.58 * (1+0.144*Mach^2)^0.65);

% flow transitions from laminar to turbulent at 500,000
if ReWinglet < ReTransition
    % if Reynolds number is less than 500k, fully laminar flow
    CfWingletSmooth = CfWingletSmoothLaminar;
else
    % Reynolds number above 500k, first part of wing in laminar flow while
    % the rest is in turbulent conditions
    CfWingletSmoothProportionLaminar = ReTransition/ReWinglet;
    CfWingletSmooth = CfWingletSmoothLaminar*CfWingletSmoothProportionLaminar + (1-CfWingletSmoothProportionLaminar)*CfWingletSmoothTurbulent;
end

% Or, with a rough surface, the following approximation from Schlichting
% CfWingletRough = (1.89 + (1.62 * log10(c/roughnessWing)))^-2.5;
CfWingletRough = 0; % NOT USED

% The greater of the two
CfWinglet = max([CfWingletSmooth CfWingletRough]);

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

R_wf = 1.0; % NO wing-fuselage interference (otherwise: Roskam page 24, taking Re_fus ~ 2-3e6 and M<0.25)
R_ls = 1.07; % Roskam page 24, based on zero sweep and M<0.25

SwetWinglet = 2*Swinglet*1.05;

FFwinglet = (1 + Lprime*(tcRatio) + 100*(tcRatio)^4);

DAwinglets = 2 * R_wf*R_ls * CfWinglet * FFwinglet * SwetWinglet;
% multiply by 2 to account for both winglets


%% Tail
% calculate laminar and turbulent skin friction coefficients
% Prandtl
CfTailSmoothLaminar = 1.328/sqrt(ReTail);
% Roskam
CfTailSmoothTurbulent = 0.455 / (log10(ReTail)^2.58 * (1+0.144*Mach^2)^0.65);

% flow transitions from laminar to turbulent at 500,000
if ReTail < ReTransition
    % if Reynolds number is less than 500k, fully laminar flow
    CfTailSmooth = CfTailSmoothLaminar;
else
    % Reynolds number above 500k, first part of wing in laminar flow while
    % the rest is in turbulent conditions
    CfTailSmoothProportionLaminar = ReTransition/ReTail;
    CfTailSmooth = CfTailSmoothLaminar*CfTailSmoothProportionLaminar + (1-CfTailSmoothProportionLaminar)*CfTailSmoothTurbulent;
end

% Or, with a rough surface, the following approximation from Schlichting
% CfTailRough = (1.89 + (1.62 * log10(c/roughnessTail)))^-2.5;
CfTailRough = 0; % NOT USED

% The greater of the two
CfTail = max([CfTailSmooth CfTailRough]);



% Form factor, Roskam (1.6 is L' for NACA 0010 with (t/c)max at 30%c
FFtail = 1 + (1.6*(tTail/cTail)) + (100*((tTail/cTail)^4));

% Total tail drag area;
DAtail = 2.05*(cTail*bTail*R_ls + cTail*hTail)*CfTail*FFtail;
% accounting for lifting surface correction, Roskam page 24


%% Fuselage

% Prandtl
CfFuseSmoothLaminar = 1.328/sqrt(ReFuse);
% Roskam
CfFuseSmoothTurbulent = 0.455 / (log10(ReFuse)^2.58 * (1+0.144*Mach^2)^0.65);

% flow transitions from laminar to turbulent at 500,000
if ReFuse < ReTransition
    % if Reynolds number is less than 500k, fully laminar flow
    CfFuseSmooth = CfFuseSmoothLaminar;
else
    % Reynolds number above 500k, first part of wing in laminar flow while
    % the rest is in turbulent conditions
    CfFuseSmoothProportionLaminar = ReTransition/ReFuse;
    CfFuseSmooth = CfFuseSmoothLaminar*CfFuseSmoothProportionLaminar + (1-CfFuseSmoothProportionLaminar)*CfFuseSmoothTurbulent;
end

% Or, with a rough surface, the following approximation from Schlichting
% CfFuseRough = (1.89 + (1.62 * log10(c/roughnessFuse)))^-2.5;
CfFuseRough = 0; % NOT USED

% The greater of the two
CfFuse = max([CfFuseSmooth CfFuseRough]);


% https://www.web-formulas.com/Math_Formulas/Geometry_Surface_of_Ellipsoid.aspx
SwetFuse = 4*pi * (((lFuse*wFuse/4)^1.6075 + (lFuse*hFuse/4)^1.6075 + (hFuse*wFuse/4)^1.6075)/3) ^ (1/1.6);
% using exponent of 1.6075 to minimize error
% dividing by 4 in each term since using 2*principal semiaxes

% Mean diameter
dFuse = 2*((wFuse*hFuse/pi)^.5);

% form factor, Roskam page 44
fuseFineness = lFuse/dFuse;
FFfuse = (1 + 60/fuseFineness^3 + 0.0025*fuseFineness);

% drag area, Roskam page 44, ignoring C_D_b_fus (small contribution)
DAfuse = R_wf * CfFuse * FFfuse * SwetFuse;


%% Landing Gear
if strcmp(lgType, 'bow')
    CDlg = 0.1;
end
Alg = hWing*(hWing/10);
DAlg = CDlg*Alg;

%% Total
DA = [DAwing DAwinglets DAfuse DAtail DApayloads DAlg];
dragResults.DAcontributors = ["Wing" "Winglets" "Fuse" "Tail" "Payloads" "Gear"];

DAtotal = sum(DA);
dragResults.DAtotal = DAtotal;
dragResults.DAshares = DA/DAtotal;

dragResults.CD0 = DAtotal/S;

if showDragBuildup % printing results
    fprintf('%8s: %.0f\n', 'Mission', mission)
    fprintf('%8s: %.3f\n', 'S', S)

    for i = 1:length(dragResults.DAcontributors)
        fprintf('%8s: %.4f -> %%%.2f of total\n', dragResults.DAcontributors(i), DA(i), dragResults.DAshares(i))
    end

    fprintf('%8s: %.4f\n', 'DA Total', DAtotal)
    fprintf('%8s: %.4f\n', 'CD0', dragResults.CD0)
    fprintf('%8s: %.3f\n', 'e', dragResults.e)

    fprintf('Press any button to continue or ctrl+break to end code...\n\n')
    pause;
end

end