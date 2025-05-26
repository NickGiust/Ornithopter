function [massFuselage] = calcMass_Fuselage(plane)
%CALCMASS_FUSELAGE Summary of this function goes here
%   Detailed explanation goes here

%% INPUTS
% epoxy
epoxy_ratio = 0.55;

% density of materials
density_carbonFiber = 5.5 / 29.494; %oz/yd2 --> kg/m2
density_fiberGlass = 3 / 29.494; %oz/yd2 --> kg/m2
density_foam = 48.06; % kg/m3 (3 lbs/ft3)
density_balsa = 160; %kg/m3
density_ply = 680; %kg/m3

% assign params from plane
wFuse = plane.wFuse;
hFuse = plane.hFuse;
lFuse = plane.lFuse;

% surface area
pEllipsoid = 1.6075;
aEllipsoid = lFuse/2;
bEllipsoid = hFuse/2;
cEllipsoid = wFuse/2;
SA_fuse = 4*pi* ( (aEllipsoid.^pEllipsoid.*bEllipsoid.^pEllipsoid + ...
    aEllipsoid.^pEllipsoid.*cEllipsoid.^pEllipsoid + ...
    bEllipsoid.^pEllipsoid.*cEllipsoid.^pEllipsoid) / 3 ).^(1/pEllipsoid);

% CF - 1 ply
m_CF = SA_fuse * density_carbonFiber;

% Wet tape joining two fuselage halves
area_wetTape = 0.25 * wFuse * 2 * lFuse;
m_wetTape = area_wetTape * density_carbonFiber;

% Fiber Glass - 1 ply
m_FG = SA_fuse * density_fiberGlass;

% Epoxy
m_epoxy_skin = (epoxy_ratio/(1-epoxy_ratio)) * (m_FG+m_CF+m_wetTape);

% Platform
area_platform = wFuse * (lFuse - 3*hFuse); % [m^2] width times (length - tailcone length)
thick_platform = 0.25 * 0.0254; % [in -> m]
m_platform = area_platform * thick_platform * density_ply; % kg

% Foam
% thickness_foam = 0.00635; %m = 1/4 inch
% 
% foam_CornerStrips_Volume = 4 * thickness_foam * lFuse * 0.2 * wFuse;
% foam_UnderWing_Volume = 2 * thickness_foam * plane.c * wFuse;
% foam_TotalVolume = foam_UnderWing_Volume + foam_CornerStrips_Volume;
% 
% m_foam = foam_TotalVolume * density_foam;

massFuselage = sum([m_CF m_wetTape m_FG m_epoxy_skin m_platform]);
end

