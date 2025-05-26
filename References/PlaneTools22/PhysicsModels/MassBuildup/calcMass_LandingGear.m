function [mLandingGear] = calcMass_LandingGear(plane)
%CALCMASS_LANDINGGEAR Summary of this function goes here
%   Detailed explanation goes here

prop2DiamInches = str2num(plane.propType2(strfind(plane.propType2,'_')+1:strfind(plane.propType2,'x')-1)); % input
prop3DiamInches = str2num(plane.propType3(strfind(plane.propType3,'_')+1:strfind(plane.propType3,'x')-1)); % input

propDiamInches = max([prop2DiamInches prop3DiamInches]);


end

