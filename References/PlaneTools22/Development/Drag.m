clear;

%% Drag Buildup Script for AF0010 Aircraft
% Jackson Markow
% 6/4/2021

% Based on method shown in Hoerner, Fluid-Dynamic Drag, Chapter 14, Sec. 2
% p 268
% * indicates subject to design change

%% PARAMETERS - ****ALL VALUES IN SI BASE UNITS****

% Airframe
b = 1.88;
ctip = .28;
croot = .28;
t = .03;
roughness = .0005;

% Components
Amotor = .0003;
Awinglet = .06;
Dnacelle = .05;
Dnacelle = 0;
propArea = .004;
gapWidth = .005;

% Conditions
Vapprox = 10;
dens = 1.2;

c = (ctip+croot)/2;
taperRatio = ctip/croot;
AR = b/c;
S = b*c;

%% Oswald efficiency
%Initial rectangular planform
f = .0524*(taperRatio^4) - 0.15*(taperRatio^3) + .1659*(taperRatio^2) - .0706*taperRatio + .0119;
eInitial = 1/(1 + (f*AR));

%Sweep
sweepEffect = .92;

%Final
e = eInitial*sweepEffect;

%% Wing surface
% Reynolds Number
Re = (dens*Vapprox*c)/1.8e-5;

% Assuming fully turbulent flow, Prandtl's power law for low Re skin
% friction
CfRe = .074/(Re^0.2);

% Or, with a rough surface, the following approximation from Schlichting
CfRough = (1.89 + (1.62 * log10(c/roughness)))^-2.5;

% The greater of the two
Cf = max([CfRe CfRough]);

% Form factor
FF = 1 + (2*(t/c)) + (60*((t/c)^4));

% Total wing drag area;
DAwing = 2.05*S*Cf*FF;

%% Other components
%Winglets
DAwinglets = 4*Awinglet*Cf;

%battery nacelle*
DAnacelle = .065*Dnacelle*c;

%Motor/mount
DAmotor = 1.15*Amotor;

%Gaps
DAgaps = .04*gapWidth*b;

%Propeller
DAprop = .09*propArea;

%% Total
DA = [DAwing DAnacelle DAmotor DAgaps DAwinglets DAprop];
DAtotal = sum(DA);
DAshares = DA/DAtotal

explode = [0 0 1 1 1 1];
DAlabels = {'Wing', 'Nacelle', 'Motor', 'Gaps', 'Winglets', 'Propeller'};
pie(DAshares, explode, DAlabels);

CD0 = DAtotal/S;