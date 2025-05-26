clc;clear; close all
fileName = "goe226-il.csv";
% read x and y data from file
    x = csvread(fileName, 9, 0, [9, 0, 41, 0]);
    y = csvread(fileName, 9, 1, [9, 1, 41, 1]);
    x = x * 0.01; % convert to m
    y = y * 0.01;

%% INPUTS
chord =  0.21; %m
span = 1.22; %m
fuselage_width = 0.13; %m

% ribs
rib_cutout = 0.8;
rib_spacing = (span/chord)*0.0117; %m
thickness_ribs = 0.003175; %m
num_ribs = floor(span/rib_spacing) + 1;
num_balsa_ribs = floor(0.5*num_ribs);
num_ply_ribs = ceil(0.5*num_ribs);

x = x * chord;
y = y * chord;

%spar
thickness_aftSpar = 0.003175; %m
thickness_frontSpar = 0.003175; %m

%spar cap
width_sparCap = 0.020; %m
num_sparCap = 8;

% epoxy
epoxy_ratio = 0.55;

density_balsa = 160; %kg/m3
density_ply = 680; %kg/m3
density_carbonFiber = 5.5 / 29.494; %oz/yd2 --> kg/yd2
density_solite = 0.7 / 29.494; %oz/yd2 --> kg/yd2

%% BUILT-UP
%% FULL AIRFOIL
% airfoil plot
airfoilPoly=polyshape(x,y);
figure
plot(airfoilPoly)
xlim([-0.2 0.6])
ylim([-0.2 0.2])
hold on

% calcs
area_fullAirfoil = polyarea(x, y);
perim_fullAirfol = perimeter(airfoilPoly);

%% CONTROL SURFACE (AFT 30%)

polyBox_controlSurface = polyshape([.7*chord chord+5 chord+5 .7*chord],[max(y) max(y) min(y) min(y)]);
plot(polyBox_controlSurface)

figure
controlSurface = intersect(airfoilPoly,polyBox_controlSurface);
plot(controlSurface)


area_ControlSurface = area(controlSurface);
perim_ControlSurface = perimeter(controlSurface);

volume_ControlSurface = area_ControlSurface*(span-fuselage_width);
m_ControlSurface = volume_ControlSurface * density_balsa;

%% RIBS (FRONT 70%)

figure
plot(airfoilPoly)
xlim([-0.2 0.6])
ylim([-0.2 0.2])
hold on
polyBox_ribs = polyshape([-1 .7*chord .7*chord -1],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
plot(polyBox_ribs)

figure
ribs = intersect(airfoilPoly,polyBox_ribs);
plot(ribs)

area_Ribs = area(ribs);
perim_Ribs = perimeter(ribs);

volume_Ply_Ribs = area_Ribs * thickness_ribs * num_ply_ribs ;
m_Ply_Ribs = density_ply * volume_Ply_Ribs;

volume_Balsa_Ribs = area_Ribs * thickness_ribs * num_balsa_ribs ;
m_Balsa_Ribs = density_balsa * volume_Balsa_Ribs;

m_Ribs = (m_Ply_Ribs + m_Balsa_Ribs) * rib_cutout;
%% FRONT SPAR

figure
plot(airfoilPoly)
xlim([-0.2 0.6])
ylim([-0.2 0.2])
hold on
polyBox_frontSpar = polyshape([.25*chord-thickness_frontSpar/2 .25*chord+thickness_frontSpar/2 .25*chord+thickness_frontSpar/2 .25*chord-thickness_frontSpar/2],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
plot(polyBox_frontSpar)

figure
frontSpar = intersect(airfoilPoly,polyBox_frontSpar);
plot(frontSpar)

area_frontSpar = area(frontSpar);
perim_frontSpar = perimeter(frontSpar);

volume_frontSpar = area_frontSpar * span;
m_frontSpar = volume_frontSpar * density_ply;
%% AFT SPAR

figure %7
plot(airfoilPoly)
xlim([-0.2 0.6])
ylim([-0.2 0.2])
hold on
polyBox_aftSpar = polyshape([.7*chord-thickness_aftSpar .7*chord .7*chord .7*chord-thickness_aftSpar],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
plot(polyBox_aftSpar)

figure %8
aftSpar = intersect(airfoilPoly,polyBox_aftSpar);
plot(aftSpar)

area_aftSpar = area(aftSpar);
perim_aftSpar = perimeter(aftSpar);

volume_aftSpar = area_aftSpar * span;
m_aftSpar = volume_aftSpar * density_ply;

%% SPAR CAPS
area_sparCap = width_sparCap * span * num_sparCap;
m_sparCap = area_sparCap * density_carbonFiber;

m_epoxy_sparCap = (epoxy_ratio/(1-epoxy_ratio)) * m_sparCap;

%% LE SKIN TO 25%

figure
plot(airfoilPoly)
xlim([-0.2 0.6])
ylim([-0.2 0.2])
hold on
polyBox_LEskin = polyshape([-1 .25*chord .25*chord -1],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
plot(polyBox_LEskin)

figure
LEskin = intersect(airfoilPoly,polyBox_LEskin);
plot(LEskin)

perim_LEskin = perimeter(LEskin);
perim_LEskin_inboard = perim_Ribs;

area_LEskin = perim_LEskin * span;
area_LEskin_inboard = perim_LEskin_inboard * fuselage_width * 2;
m_LEskin = (area_LEskin + area_LEskin_inboard) * density_carbonFiber;

m_epoxy_LEskin = (epoxy_ratio/(1-epoxy_ratio)) * m_LEskin;

%% SOLITE 25% -->
figure
plot(airfoilPoly)
xlim([-0.2 0.6])
ylim([-0.2 0.2])
hold on
polyBox_TEskin = polyshape([.25*chord chord+5 chord+5 .25*chord],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
plot(polyBox_TEskin)

figure
TEskin = intersect(airfoilPoly,polyBox_TEskin);
plot(TEskin)

perim_TEskin = perimeter(TEskin);

area_TEskin = perim_TEskin * (span - 2*fuselage_width);
m_TEskin = area_TEskin * density_solite;

m_builtupWing = sum([m_ControlSurface m_Ribs m_frontSpar m_aftSpar m_sparCap m_epoxy_sparCap m_LEskin m_epoxy_LEskin m_TEskin]);



m_builtupWing
m_actual = .484 - 0.044
delta = m_builtupWing - m_actual
