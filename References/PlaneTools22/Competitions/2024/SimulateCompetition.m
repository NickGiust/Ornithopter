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

% TOFL setup
TOFLLimit = 20 * 0.3048; % [ft -> m] from the rules (20 ft)

%% MISSION 2
% Build up weight and drag
dragResults2 = dragBuildup(plane, 2);
plane.CD02 = dragResults2.CD0;
plane.e = dragResults2.e;
plane.dragResults2 = dragResults2;

% Set mission parameters
mission.missionNo = 2;
mission.lapLimit = 3;
mission.timeLimit = 300; % 5 minutes
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.climbCL = 0.8; % JTS - update with # from David
mission.turnCL = 0.7*plane.CLmax - 0.1; % non-max CL to not account for flaps - 0.1 (to account for push down from tail)

crzThrottleList2 = 1; % sets max value for the throttle to be set at
missionResults2Array = [];
timeArray = zeros(1, length(crzThrottleList2)); % preallocating for speed

for i = 1:length(crzThrottleList2)
    mission.crzThrottle = crzThrottleList2(i);
    mission.turnThrottle = mission.crzThrottle;
    TOthrottleConvergence = false;
    missionResults2 = SimulateMission(plane, environment, mission, TOFLLimit, TOthrottleConvergence);
    missionResults2Array = [missionResults2Array missionResults2]; % JTS - try to preallocate
    timeArray(i) = missionResults2.missionTime(end);
end

[~, minTimeIndex] = min(timeArray);

missionResults2 = missionResults2Array(minTimeIndex);
missionResults2.optimizedThrottleSetting = crzThrottleList2(minTimeIndex);


%% MISSION 3
% Build up weight and drag
plane.m3 = plane.m2 - plane.mPayload2 + plane.mPayload3;
dragResults3 = dragBuildup(plane, 3);
plane.CD03 = dragResults3.CD0;
plane.dragResults3 = dragResults3;

% Set mission parameters
mission.missionNo = 3;
mission.lapLimit = 40;
mission.timeLimit = 300; % 5 minutes
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.turnThrottle = 0.8;
mission.climbCL = 0.8; % JTS - update with # from David
mission.turnCL = 0.8*plane.CLmax; % JTS - update with # from David

mission.crzThrottle = 0.8;

if missionResults2.TOFLbool == "FAILED" || missionResults2.missionTime <= 0
    % Don't run mission 3 if mission 2 failed
    missionResults3.missionNo = 3;
    missionResults3.TOFLbool = "NOT ATTEMPTED";
    missionResults3.TOFL = 0;
    missionResults3.VTO = 0;
    missionResults3.nLaps = 0;
    missionResults3.missionTime = 0;
    missionResults3.minTOthrottle = 1.1;
    missionResults3.optimizedThrottleSetting = 0;
else
    % Simulate mission
    crzThrottleList3 = 1;
    missionResults3Array = [];
    nLapsArray = zeros(1, length(crzThrottleList3)); % preallocating for speed
    
    for i = 1:length(crzThrottleList3)
        mission.crzThrottle = crzThrottleList3(i);
        mission.turnThrottle = mission.crzThrottle;
        TOthrottleConvergence = false;
        missionResults3 = SimulateMission(plane, environment, mission, TOFLLimit, TOthrottleConvergence);
        missionResults3Array = [missionResults3Array missionResults3]; % JTS - try to preallocate
        nLapsArray(i) = missionResults3.nLaps;
    end
    
    [~, maxLapsIndex] = max(nLapsArray);
    
    missionResults3 = missionResults3Array(maxLapsIndex);
    missionResults3.optimizedThrottleSetting = crzThrottleList3(maxLapsIndex);
end

end

