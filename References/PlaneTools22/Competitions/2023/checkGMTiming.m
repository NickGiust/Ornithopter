function [trueGMloadingMultiplier] = checkGMTiming(plane)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% times (seconds)
takingPartsOut = 20;
mainGear = 15;
tailBoom = 10;
package = 15;
vstab = 10;
hstab = 15;
tailWiring = 15;
centerFairings = 20;
leftWing = 20;
leftWingWiring = 10;
rightWing = leftWing;
rightWingWiring = leftWingWiring;
fixture = 90;
straps = 30;
hoistLifting = 60;
hold = 0;

setupTime = sum([takingPartsOut mainGear tailBoom package vstab hstab tailWiring centerFairings leftWing leftWingWiring rightWing rightWingWiring fixture straps hoistLifting hold]);

remainingTime = (10 * 60) - setupTime; % [seconds]; 10 min window

timePerMass = 10/(60/2.205); % 10s to load 60 lb (27.22 kg) weight

loadableMass = remainingTime/timePerMass;

maxGMloadingMultiplier = loadableMass/plane.m2;

trueGMloadingMultiplier = min([maxGMloadingMultiplier plane.GMloadingMultiplier]);

end