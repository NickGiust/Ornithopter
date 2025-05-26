function [missionResults2, missionResults3, missionResultsGM, plane] = TS_SimulateCompetition(plane, environment, planeToolsDirectory)

% global planeToolsDirectory;

if ~plane.planeFits_boolean || ~plane.antennaFits_boolean
    missionResults2.score = 0;
    missionResults2.TOFLbool = "NOT ATTEMPTED";
    missionResults3.score = 0;
    missionResults3.TOFLbool = "NOT ATTEMPTED";
    missionResultsGM.score = 0;
    %fprintf('BOX DID NOT FIT')
    return
end


%% Set up thrust model

global T CTf CPf
CTf = 1; % multiplier on thrust coefficient (empirical/theoretical --> experimental)
CPf = 1; % multiplier on thrust coefficient (empirical/theoretical --> experimental)
if ~isa(T,'function_handle') 
    ThrustModel; 
end

%% MASS BUILD UP LOOP

% Build up weight and drag
plane.e = dragBuildup(plane, 2).e;
plane.seed_m2 = 6; %kg, aircraft weight guess
currentRatio = 0;
alpha = 0.5;
num_MassBuildups = 0;

while abs(currentRatio - plane.GMloadingMultiplier) > (0.02*plane.GMloadingMultiplier)
    [M2massResults, plane] = massBuildup(plane, 2);
    currentRatio = (plane.seed_m2*plane.GMloadingMultiplier) / M2massResults.m;
    seed_m2 = plane.seed_m2 - alpha * plane.seed_m2 * (currentRatio - plane.GMloadingMultiplier) / plane.GMloadingMultiplier;
    plane.seed_m2 = seed_m2;
    num_MassBuildups = num_MassBuildups + 1;
    if num_MassBuildups >= 40 || plane.validTubes == false
        missionResults2.score = 0;
        missionResults3.score = 0;
        missionResultsGM.score = 0;
        missionResults2.TOFLbool = "NOT ATTEMPTED";
        fprintf('\nMASS COULD NOT CONVERGE\n\n')
        return
    end
end
fprintf("Mass convergence loops: %g\n", num_MassBuildups)
plane.GMloadingMultiplier = currentRatio;

plane.m2 = M2massResults.m;
plane.m2massResults = M2massResults;
plane.trueGMloadingMultiplier = checkGMTiming(plane);
missionResultsGM.score = plane.trueGMloadingMultiplier;

%% MISSION 2
% Build up weight and drag
dragResults2 = dragBuildup(plane, 2);
plane.CD02 = dragResults2.CD0;
plane.dragResults2 = dragResults2;

plane.CD02 = dragBuildup(plane, 2).CD0;

% Set mission parameters
mission.missionNo = 2;
mission.lapLimit = 40;
mission.timeLimit = 600; % 10 minutes
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.crzThrottle = 0.53;
mission.turnThrottle = 0.53;
mission.climbCL = 1.0;
mission.turnCL = 0.85*plane.CLmax;

TOFLLimit = 60*0.3048;

missionResults2 = SimulateMission(plane, environment, mission, TOFLLimit);

% PlotPerformance(mission, plane, missionResults2);
% CheckLimits(mission, plane, missionResults2);

%% MISSION 3
% Build up weight and drag
plane.m3 = plane.m2 - (plane.payload2Fraction * plane.m2) + (2 * plane.lengthPayload3 * 0.238106231 + 0.1);
dragResults3 = dragBuildup(plane, 3);
plane.CD03 = dragResults3.CD0;
plane.dragResults3 = dragResults3;


% Set mission parameters
mission.missionNo = 3;
mission.lapLimit = 3;
mission.timeLimit = 300; % 5 minutes
mission.TOthrottle = 0.8;
mission.climbThrottle = 0.8;
mission.crzThrottle = 0.9;
mission.turnThrottle = 0.9;
mission.climbCL = 0.6;
mission.turnCL = 0.8*plane.CLmax;

if missionResults2.TOFLbool == "FAILED"
    missionResults3.score = 0;
    missionResults3.missionNo = 3;
    missionResults3.TOFLbool = "NOT ATTEMPTED";
    missionResults3.TOFL = 0;
    missionResults3.nLaps = 0;
    missionResults3.missionTime = 0;
else
    % Simulate mission
    missionResults3 = SimulateMission(plane, environment, mission, TOFLLimit);
end

% PlotPerformance(mission, plane, missionResults3);
% CheckLimits(mission, plane, missionResults3);





end

