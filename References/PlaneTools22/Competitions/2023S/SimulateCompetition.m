function [missionResults2, missionResults3, plane] = SimulateCompetition(plane, environment)

%% Set up thrust model

global T CTf CPf
CTf = 1; % multiplier on thrust coefficient (empirical/theoretical --> experimental)
CPf = 1; % multiplier on thrust coefficient (empirical/theoretical --> experimental)
if ~isa(T,'function_handle') 
    ThrustModel; 
end

% Build up weight and drag
[M2massResults, plane] = massBuildup(plane, 2);
plane.m2 = M2massResults.m;

%% MISSION 2
% Build up weight and drag
dragResults2 = dragBuildup(plane, 2);
plane.CD02 = dragResults2.CD0;
plane.e = dragResults2.e;
plane.dragResults2 = dragResults2;

% Set mission parameters
mission.missionNo = 2;
mission.lapLimit = 40;
mission.timeLimit = 600; % 10 minutes
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.climbCL = 1.0;
mission.turnCL = 0.85*plane.CLmax;

TOFLLimit = 25 * 0.3048; % [ft -> m]
crzThrottleList = 0.3:0.05:0.6;
missionResults2Array = [];
nLapsArray = zeros(1, length(crzThrottleList)); % preallocating for speed

for i = 1:length(crzThrottleList)
    mission.crzThrottle = crzThrottleList(i);
    mission.turnThrottle = mission.crzThrottle;
    TOthrottleConvergence = true;
    missionResults2 = SimulateMission(plane, environment, mission, TOFLLimit, TOthrottleConvergence);
    missionResults2Array = [missionResults2Array missionResults2]; % JTS - try to preallocate
    nLapsArray(i) = missionResults2.nLaps;
end

[~, maxLapsIndex] = max(nLapsArray);

missionResults2 = missionResults2Array(maxLapsIndex);
missionResults2.optimizedThrottleSetting = crzThrottleList(maxLapsIndex);


%% MISSION 3
% Build up weight and drag
plane.m3 = plane.m2 - plane.radarMass + plane.nPayloads3 * plane.m3DropsondeMass;
dragResults3 = dragBuildup(plane, 3);
plane.CD03 = dragResults3.CD0;
plane.dragResults3 = dragResults3;

% Set mission parameters
mission.missionNo = 3;
mission.lapLimit = 3;
mission.timeLimit = 300; % 5 minutes
mission.TOthrottle = 0.95;
mission.climbThrottle = 0.95;
mission.crzThrottle = 0.8;
mission.turnThrottle = 0.8;
mission.climbCL = 0.6; % JTS - update with # from David
mission.turnCL = 0.8*plane.CLmax; % JTS - update with # from David

if missionResults2.TOFLbool == "FAILED"
    % Don't run mission 3 if mission 2 failed
    missionResults3.missionNo = 3;
    missionResults3.TOFLbool = "NOT ATTEMPTED";
    missionResults3.TOFL = 0;
    missionResults3.nLaps = 0;
    missionResults3.missionTime = 0;
    missionResults3.minTOthrottle = 1.1;
else
    % Simulate mission
    TOthrottleConvergence = true;
    missionResults3 = SimulateMission(plane, environment, mission, TOFLLimit,TOthrottleConvergence);
end

end

