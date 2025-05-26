function [massBuiltupWing, plane] = calcMass_BuiltupWing(plane)
%CALCMASS_BUILTUPWING Summary of this function goes here
%   Detailed explanation goes here

%% INPUTS

% setting params from plane
chord = plane.c;
span = plane.b;
fuselage_width = plane.wFuse;
control_frac = plane.wingControlFraction;

% reading in airfoil x, y coordinates
airfoilFileName = plane.airfoilCoordinatesFile; % need to have this as input!
    % read file (specifically the appropriate lines only)
    airfoilDataOpts = detectImportOptions(airfoilFileName);
    airfoilDataOpts.DataLines = [10, inf]; % NOTE: ASSUMES THAT THE AIRFOIL DATA STARTS ON LINE 10
    airfoilData = readmatrix(airfoilFileName, airfoilDataOpts);
    % convert from c = 100 -> c = 1
    airfoilData = airfoilData ./ 100;
    % adjust to chord length
    airfoilData = airfoilData .* chord;
    % separate into x and y
    x = airfoilData(:,1);
    y = airfoilData(:,2);

% rib assumptions
rib_cutout = 0.8; % percent of wood left after mass-reduction cuts made
thickness_ribs = 0.003175; %m = 1/8 inch
num_ribs = 10 * 2; % per side -> total
rib_spacing = (span - plane.wFuse - 2*thickness_ribs)/(2*num_ribs - 2); %m, space between ribs on one side
% num_ply_ribs = ceil(0.6*num_ribs);
num_ply_ribs = 3 * 2; % per side -> total, 1 each for flaps, ailerons, base
% num_balsa_ribs = floor(0.4*num_ribs);
num_balsa_ribs = num_ribs - num_ply_ribs;
epoxy_glue = 0.01; %kg, weight of epoxy used per rib

% spar assumptions
thickness_mainSpar = 0.003175; %m = 1/8 in
thickness_aftSpar = 0.003175; %m = 1/8 inch

% spar cap (INTEGRATE SPAR SIZER HERE)
% [tubePlyArea, shearWebVolume, plane] = tubeSizer(plane);

% epoxy
epoxy_ratio = 0.55;

% density of materials
density_balsa = 160; %kg/m3
density_ply = 680; %kg/m3
density_foam = 3 * 16.018; % converting 3 lb/ft^3 foam to kg/m^3
density_carbonFiber = 5.5 / 29.494; %oz/yd2 --> kg/m2
density_uniCarbonFiber = 4.7 / 29.494; %oz/yd2 --> kg/m2 from https://www.cstsales.com/uni_carbon_fabric-ss2.html
density_solite = 0.7 / 29.494; %oz/yd2 --> kg/m2

%% FULL AIRFOIL
airfoilPoly=polyshape(x,y);
plane.airfoilPerimeter = perimeter(airfoilPoly);

%% CONTROL SURFACE (AFT 30%)
polyBox_controlSurface = polyshape([(1-control_frac)*chord chord+5 chord+5 (1-control_frac)*chord],[max(y) max(y) min(y) min(y)]);
controlSurface = intersect(airfoilPoly,polyBox_controlSurface);
area_ControlSurface = area(controlSurface);
volume_ControlSurface = area_ControlSurface*(span-fuselage_width);

if strcmpi(plane.controlSurfMat, "balsa")
    m_ControlSurface = volume_ControlSurface * density_balsa;
elseif strcmpi(plane.controlSurfMat, "foam")
    m_ControlSurface = volume_ControlSurface * density_foam;
end

%% RIBS (FRONT 70%)
polyBox_ribs = polyshape([-1 (1-control_frac)*chord (1-control_frac)*chord -1],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
ribs = intersect(airfoilPoly,polyBox_ribs);
area_Ribs = area(ribs);
perim_Ribs = perimeter(ribs);

volume_Ply_Ribs = area_Ribs * thickness_ribs * num_ply_ribs ;
m_Ply_Ribs = density_ply * volume_Ply_Ribs;

volume_Balsa_Ribs = area_Ribs * thickness_ribs * num_balsa_ribs ;
m_Balsa_Ribs = density_balsa * volume_Balsa_Ribs;

m_Ribs = (m_Ply_Ribs + m_Balsa_Ribs) * rib_cutout;
m_epoxy_ribs = epoxy_glue * num_ribs;

%% MAIN SPAR

polyBox_mainSpar = polyshape([.25*chord-thickness_mainSpar/2 .25*chord+thickness_mainSpar/2 .25*chord+thickness_mainSpar/2 .25*chord-thickness_mainSpar/2],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
mainSpar = intersect(airfoilPoly,polyBox_mainSpar);
area_mainSpar = area(mainSpar);

volume_mainSpar = area_mainSpar * span;
m_mainSpar = volume_mainSpar * density_ply;

%% AFT SPAR
polyBox_aftSpar = polyshape([(1-control_frac)*chord-thickness_aftSpar (1-control_frac)*chord (1-control_frac)*chord (1-control_frac)*chord-thickness_aftSpar],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
aftSpar = intersect(airfoilPoly,polyBox_aftSpar);
area_aftSpar = area(aftSpar);

volume_aftSpar = area_aftSpar * span;
m_aftSpar = volume_aftSpar * density_ply;

%% SPAR CAPS
%m_sparCap = tubePlyArea * density_uniCarbonFiber;

%m_epoxy_sparCap = (epoxy_ratio/(1-epoxy_ratio)) * m_sparCap;

%% LE SKIN TO dboxFraction (as percent)
dboxFraction = 0.30;
LEskin_inboard_width = 2.5 * fuselage_width;
polyBox_LEskin = polyshape([-1 dboxFraction*chord dboxFraction*chord -1],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
LEskin = intersect(airfoilPoly,polyBox_LEskin);
perim_LEskin = perimeter(LEskin);
perim_LEskin_inboard = perim_Ribs;

area_LEskin = perim_LEskin * (span-LEskin_inboard_width);
area_LEskin_inboard = perim_LEskin_inboard * LEskin_inboard_width;
m_LEskin = (area_LEskin + area_LEskin_inboard) * density_carbonFiber;

m_epoxy_LEskin = (epoxy_ratio/(1-epoxy_ratio)) * m_LEskin;

%% SOLITE 25% onward -->
polyBox_TEskin = polyshape([dboxFraction*chord chord+5 chord+5 dboxFraction*chord],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
TEskin = intersect(airfoilPoly,polyBox_TEskin);
perim_TEskin = perimeter(TEskin);

area_TEskin = perim_TEskin * (span - LEskin_inboard_width);

if strcmpi(plane.controlSurfMat, "balsa")
    m_TEskin = area_TEskin * density_solite;
elseif strcmpi(plane.controlSurfMat, "foam")
    m_TEskin = area_TEskin * density_fiberGlass;
end

%% ADD IT ALL UP
%wingMassComponents = [m_ControlSurface m_Ribs m_mainSpar m_aftSpar m_sparCap m_epoxy_sparCap m_LEskin m_epoxy_LEskin m_TEskin];
wingMassComponents = [m_ControlSurface m_Ribs m_epoxy_ribs m_mainSpar m_aftSpar m_LEskin m_epoxy_LEskin m_TEskin];

massBuiltupWing = sum(wingMassComponents);

end

