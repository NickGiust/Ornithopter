function [mESC] = calcMass_ESC(plane, ESC_continuousCurrent)
%CALCMASS_ESC Summary of this function goes here
%   Detailed explanation goes here

% pheonix edge lite series (castle creations)

currentMargin = 1.4; % how far above current rating we willin to risk it
ESC_allowableCurrent = ESC_continuousCurrent/currentMargin;

ESC_rating = [50,75,100,130,200];
ESC_mass = [56,81,91.5,115,187]/1000;

index = find(ESC_rating > ESC_allowableCurrent, 1);
if isempty(index)
    index = length(ESC_rating);
end

mESC = ESC_mass(index) * plane.nMotors;
end

