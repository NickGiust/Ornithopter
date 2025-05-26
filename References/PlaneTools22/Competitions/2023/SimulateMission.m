function missionResults = SimulateMission(plane, environment, mission, TOFLLimit, TOthrottleConvergence)
%% Mission parameters
missionNo = mission.missionNo;
lapLimit = mission.lapLimit;
timeLimit = mission.timeLimit;
TOthrottle = mission.TOthrottle;
climbThrottle = mission.climbThrottle;
crzThrottle = mission.crzThrottle;
turnThrottle = mission.turnThrottle;
climbCL = mission.climbCL;
turnCL = mission.turnCL;

segment.coursePoints.segments = [];
segment.coursePoints.startTimes = [];
nLaps = 0;
segment.performance = [];
missionPoints.lapNos = 1;
missionPoints.lapStartTimes = 0;
timeAtLapEnd = 0;
decelDistance = 0;
warning('off','all');

%% THE MISSION
% Takeoff

if TOthrottleConvergence
    TOthrottleList = 0.7:0.1:1;
else
    TOthrottleList = 1;
end
for i = 1:length(TOthrottleList)
    TOthrottle = TOthrottleList(i);
    segment = FlightModel(plane,environment,missionNo,TOthrottle,0,'takeoff',0,0,segment.coursePoints,segment.performance);
    TOFL = segment.TOFL;
    VTO = segment.VTO;
    Vstall = segment.Vstall;
    missionResults.TOFL = TOFL;
    missionResults.VTO = VTO;
    missionResults.Vstall = Vstall;
    
    if TOFL >= 0 && TOFL < TOFLLimit
        missionResults.TOFLbool = "PASSED";
        minTOthrottle = TOthrottleList(i);
        break
    else
        segment.coursePoints.segments = [];
        segment.coursePoints.startTimes = [];
        segment.performance = [];
    end

end

if TOFL > TOFLLimit || TOFL <= 0
    missionResults.TOFLbool = "FAILED";
    missionResults.score = 0;
    missionResults.missionNo = missionNo;
    missionResults.nLaps = 0;
    missionResults.missionTime = 0;
    missionResults.minTOthrottle = 1.1;

%     segment.performance
%     segment.performance.t
%     segment.performance.V
%     segment.performance.s

    return
end
% segment.performance
% segment.performance.t
% segment.performance.V
% segment.performance.s

% Climb
segment = FlightModel(plane,environment,missionNo,climbThrottle,climbCL,'climb',0,0,segment.coursePoints,segment.performance);

usableBatteryFraction = 0.85; % Must be 0-0.99

%% Cruise loop
while nLaps < lapLimit && timeAtLapEnd < timeLimit && segment.performance.s(end) > (1-usableBatteryFraction)
    
    batteryChargeAtLapStart = segment.performance.s(end);

    missionPoints.lapNos = [missionPoints.lapNos; nLaps+1];
    missionPoints.lapStartTimes = [missionPoints.lapStartTimes; segment.performance.t(end)];

    % Half straight upwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',1, ...
        segment.additionalDistance,segment.coursePoints,segment.performance);
  
    % 180
    segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn180',0,0,segment.coursePoints,segment.performance);

    % Half straight downwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',0, ...
        segment.additionalDistance,segment.coursePoints,segment.performance);

    % 360
    segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn360',0,0,segment.coursePoints,segment.performance);
        
    % Half straight downwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',0, ...
        segment.additionalDistance - min([0 153+decelDistance]),segment.coursePoints,segment.performance);

    % 180
    segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn180',0,0,segment.coursePoints,segment.performance);

    % Half straight upwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',1, ...
        segment.additionalDistance + decelDistance,segment.coursePoints,segment.performance);

    timeAtLastLapEnd = timeAtLapEnd;
    timeAtLapEnd = segment.performance.t(end);
    nLaps = nLaps + 1;
end

% if segment.performance.s(end) < (1-usableBatteryFraction)
%     nLaps = nLaps - 1;
% end

if timeAtLapEnd < timeLimit
    missionTime = timeAtLapEnd;
else
    missionTime = timeAtLastLapEnd;
    nLaps = nLaps - 1;
end

%% OUTPUT TO PLANETOOLS


if missionNo == 2
    missionResults.score = plane.mPayload2 * nLaps;
elseif missionNo == 3
    missionResults.score = plane.lengthPayload3 / missionTime;
end


missionResults.missionNo = missionNo;
missionResults.TOFL = TOFL;
missionResults.minTOthrottle = minTOthrottle;
missionResults.VTO = VTO;
missionResults.nLaps = nLaps;
missionResults.missionTime = missionTime;
missionResults.missionPoints = missionPoints;
missionResults.coursePoints = segment.coursePoints;
missionResults.performance = segment.performance;
missionResults.batteryUsage = 1 - batteryChargeAtLapStart;



end


