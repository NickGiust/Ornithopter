function [massFoamcoreWing, plane] = calcMass_FoamcoreWing(plane)
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
thickness_ribs = 0.003175; % [inch -> m], 1/8"
if rem(plane.nMotors, 2) == 0
    num_ribs = 2 + plane.nMotors * 4; % if even number of motors, 2 ribs for connecting to fuse + 2 for every motor mount
else
    num_ribs = 2 + (plane.nMotors - 1) * 4; % if odd number of motors, 2 ribs for connecting to fuse + 2 for every motor ON WING
end
epoxy_glue = 0.01; %kg, weight of epoxy used per rib to glue into foamcore wing

% spar assumptions
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
density_fiberGlass = 3 / 29.494; %oz/yd2 --> kg/m2

%% FULL AIRFOIL
airfoilPoly=polyshape(x,y);
plane.airfoilPerimeter = perimeter(airfoilPoly);

%% CONTROL SURFACE (AFT portion)
polyBox_controlSurface = polyshape([(1-control_frac)*chord chord+5 chord+5 (1-control_frac)*chord],[max(y) max(y) min(y) min(y)]);
controlSurface = intersect(airfoilPoly,polyBox_controlSurface);
area_ControlSurface = area(controlSurface);
volume_ControlSurface = area_ControlSurface*(span-fuselage_width);

if strcmpi(plane.controlSurfMat, "balsa")
    m_ControlSurface = volume_ControlSurface * density_balsa;
elseif strcmpi(plane.controlSurfMat, "foam")
    m_ControlSurface = volume_ControlSurface * density_foam;
end

%% FOAMCORE (FRONT PORTION OF WING)
polyBox_core = polyshape([-1 (1-control_frac)*chord (1-control_frac)*chord -1],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
core = intersect(airfoilPoly,polyBox_core);
area_core = area(core);
perim_core = perimeter(core);

volume_core = area_core * (span - num_ribs);
m_core = density_foam * volume_core;

%% RIBS
volume_ribs = area_core * thickness_ribs * num_ribs;

m_ribs = density_ply * volume_ribs; % assumes all ribs are plywood (mostly all load bearing)
m_epoxy_ribs = epoxy_glue * num_ribs;

%% FIBERGLASS SKIN
area_skin = perim_core * span;

m_LEskin = area_skin * density_fiberGlass;
m_epoxy_LEskin = (epoxy_ratio/(1-epoxy_ratio)) * m_LEskin;

%% AFT SPAR
polyBox_aftSpar = polyshape([(1-control_frac)*chord-thickness_aftSpar (1-control_frac)*chord (1-control_frac)*chord (1-control_frac)*chord-thickness_aftSpar],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
aftSpar = intersect(airfoilPoly,polyBox_aftSpar);
area_aftSpar = area(aftSpar);

volume_aftSpar = area_aftSpar * span;
m_aftSpar = volume_aftSpar * density_ply;

%% SPAR CAPS (assumes 2" wide sparcaps made of uni on top and bottom, 1 ply each)
polyBox_sparCap = polyshape([0.25*chord-0.0254 0.25*chord+0.0254 0.25*chord+0.0254 0.25*chord-0.0254],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
sparCap = intersect(airfoilPoly,polyBox_sparCap);
perim_sparCap = perimeter(sparCap);

m_sparCap = perim_sparCap * density_uniCarbonFiber;
m_epoxy_sparCap = (epoxy_ratio/(1-epoxy_ratio)) * m_sparCap;

%% LE SKIN TO dboxFraction (as percent)
% dboxFraction = 0.30;
% LEskin_inboard_width = 2.5 * fuselage_width;
% polyBox_LEskin = polyshape([-1 dboxFraction*chord dboxFraction*chord -1],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
% LEskin = intersect(airfoilPoly,polyBox_LEskin);
% perim_LEskin = perimeter(LEskin);
% perim_LEskin_inboard = perim_Ribs;
% 
% area_LEskin = perim_LEskin * (span-LEskin_inboard_width);
% area_LEskin_inboard = perim_LEskin_inboard * LEskin_inboard_width;
% m_LEskin = (area_LEskin + area_LEskin_inboard) * density_carbonFiber;
% 
% m_epoxy_LEskin = (epoxy_ratio/(1-epoxy_ratio)) * m_LEskin;

%% FIBERGLASS CONTROL SURFACES
polyBox_TEskin = polyshape([(1-control_frac)*chord chord+5 chord+5 (1-control_frac)*chord],[max(y)+5 max(y)+5 min(y)-5 min(y)-5]);
TEskin = intersect(airfoilPoly,polyBox_TEskin);
perim_TEskin = perimeter(TEskin);

area_TEskin = perim_TEskin * (span-fuselage_width);

if strcmpi(plane.controlSurfMat, "balsa")
    m_TEskin = area_TEskin * density_solite;
elseif strcmpi(plane.controlSurfMat, "foam")
    m_TEskin = area_TEskin * density_fiberGlass;
end

%% ADD IT ALL UP
wingMassComponents = [m_ControlSurface m_core m_ribs m_epoxy_ribs m_LEskin m_epoxy_LEskin m_aftSpar m_sparCap m_epoxy_sparCap m_TEskin];

massFoamcoreWing = sum(wingMassComponents);

end

