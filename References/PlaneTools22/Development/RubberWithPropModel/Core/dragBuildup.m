function dragResults = dragBuildup(plane, mission)
%% PARAMETERS
% Wing
b = plane.b;
ctip = plane.cTip;
croot = plane.cRoot;
t = plane.wingThickness;

load('material.mat', plane.wingMat);
wingMaterial = eval(plane.wingMat);
roughnessWing = wingMaterial.roughness;

DAwinglets = plane.DAwinglets;

% Tail
cTail = plane.cTail;
bTail = plane.bTail;
hTail = plane.hTail;
tTail = cTail*.1;
Stail = cTail*(bTail+hTail);

load('material.mat', plane.tailMat);
tailMaterial = eval(plane.tailMat);
roughnessTail = tailMaterial.roughness;

% Fuselage
lFuse = plane.lFuse;
wFuse = plane.wFuse;
hFuse = plane.hFuse;

load('material.mat', plane.fuseMat);
fuselageMaterial = eval(plane.fuseMat);
roughnessFuse = fuselageMaterial.roughness;

% Payloads
if mission == 2
    DApayloads = plane.nPayloads2*plane.DApayload2;
elseif mission == 3
    DApayloads = plane.nPayloads3*plane.DApayload3;
end

% Landing gear
lgType = plane.lgType;
hWing = plane.hWing;

% Air
Vapprox = 10;
dens = 1.2;

% Basic equations
c = (ctip+croot)/2;
taperRatio = ctip/croot;
AR = b/c;
S = b*c;

%% Reynolds Numbers
% Wing
ReWing = (dens*Vapprox*c)/1.8e-5;

% Fuselage
ReFuse = (dens*Vapprox*lFuse)/1.8e-5;

% Tail 
ReTail = (dens*Vapprox*cTail)/1.8e-5;

%% Oswald efficiency
% Initial rectangular planform
f = .0524*(taperRatio^4) - 0.15*(taperRatio^3) + .1659*(taperRatio^2) - .0706*taperRatio + .0119;
eInitial = 1/(1 + (f*AR));

% Sweep
sweepEffect = .92;

% Final
dragResults.e = eInitial*sweepEffect;

%% Wing
% Assuming fully turbulent flow, Prandtl's power law for low Re skin
% friction
CfWingSmooth = .074/(ReWing^0.2);

% Or, with a rough surface, the following approximation from Schlichting
CfWingRough = (1.89 + (1.62 * log10(c/roughnessWing)))^-2.5;

% The greater of the two
CfWing = max([CfWingSmooth CfWingRough]);

% Form factor
FFwing = 1 + (2*(t/c)) + (60*((t/c)^4));

% Total wing drag area;
DAwing = 2.05*S*CfWing*FFwing;

%% Tail
% Assuming fully turbulent flow, Prandtl's power law for low Re skin
% friction
CfTailSmooth = .074/(ReTail^0.2);

% Or, with a rough surface, the following approximation from Schlichting
CfTailRough = (1.89 + (1.62 * log10(c/roughnessTail)))^-2.5;

% The greater of the two
CfTail = max([CfTailSmooth CfTailRough]);

% Form factor
FFtail = 1 + (2*(tTail/cTail)) + (60*((tTail/cTail)^4));

% Total tail drag area;
DAtail = 2.05*Stail*CfTail*FFtail;

%% Fuselage
% Assuming fully turbulent flow, Prandtl's power law for low Re skin
% friction
CfFuseSmooth = .074/(ReFuse^0.2);

% Or, with a rough surface, the following approximation from Schlichting
CfFuseRough = (1.89 + (1.62 * log10(c/roughnessFuse)))^-2.5;

% The greater of the two
CfFuse = max([CfFuseSmooth CfFuseRough]);

% Mean diameter
dFuse = 2*((wFuse*hFuse/pi)^.5);

% Form factor
FRfuse = lFuse/dFuse;
FFfuse = 1 + (1.5/(FRfuse^1.5)) + (7/(FRfuse^3));

% Wetted area
Sfuse = wFuse*hFuse*4*(FRfuse - 1.3);

% Total fuselage drag area;
DAfuse = 2.05*Sfuse*CfFuse*FFfuse;

%% Landing Gear
if strcmp(lgType, 'strut')
    CDlg = .3;
end
Alg = hWing*(hWing/10);
DAlg = CDlg*Alg;

%% Total
DA = [DAwing DAwinglets DAfuse DAtail DApayloads DAlg];
DAtotal = sum(DA);
dragResults.DAshares = DA/DAtotal;

dragResults.CD0 = DAtotal/S;
end