function missionResults = SimulateMission(plane, environment, mission)
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
segment = FlightModel(plane,environment,missionNo,TOthrottle,0,'takeoff',0,0,segment.coursePoints,segment.performance);
TOFL = segment.TOFL;
VTO = segment.VTO;
  
% Climb
segment = FlightModel(plane,environment,missionNo,climbThrottle,climbCL,'climb',0,0,segment.coursePoints,segment.performance);

%% Cruise loop
while nLaps < lapLimit && timeAtLapEnd < timeLimit && segment.performance.s(end) > .01
    
    if nLaps > 0
        missionPoints.lapNos = [missionPoints.lapNos; nLaps+1];
        missionPoints.lapStartTimes = [missionPoints.lapStartTimes; segment.performance.t(end)];
    end

    if missionNo == 3 && nLaps > 0
        % Takeoff
        segment = FlightModel(plane,environment,missionNo,TOthrottle,0,'takeoff',0,0,segment.coursePoints,segment.performance);
        TOFL = segment.TOFL;
        VTO = segment.VTO;

        % Climb
        segment = FlightModel(plane,environment,missionNo,climbThrottle,climbCL,'climb',0,0,segment.coursePoints,segment.performance);
    end
    
    if nLaps > 0
    
        % Half straight upwind
        segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',1, ...
            segment.additionalDistance,segment.coursePoints,segment.performance);
        
        if nLaps < plane.nPayloads3 && missionNo == 3
            plane.m3 = plane.m3 - plane.mPayload3;
        end
    end

    % 180
    segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn180',0,0,segment.coursePoints,segment.performance);

    % Half straight downwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',0, ...
        segment.additionalDistance,segment.coursePoints,segment.performance);

    %{
    % 360
    segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn360',0,0,segment.coursePoints,segment.performance);
    %}
        
    % Half straight downwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',0, ...
        segment.additionalDistance - min([0 153+decelDistance]),segment.coursePoints,segment.performance);

    % 180
    segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn180',0,0,segment.coursePoints,segment.performance);

    % Half straight upwind
    segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',1, ...
        segment.additionalDistance + decelDistance,segment.coursePoints,segment.performance);
    
    if missionNo == 3
        % Deceleration
        segment = FlightModel(plane,environment,missionNo,0,0,'deceleration',1, ...
            segment.additionalDistance,segment.coursePoints,segment.performance);
        decelDistance = segment.additionalDistance;

        % Stopping
        segment = FlightModel(plane,environment,missionNo,0,0,'stopping',1, ...
            segment.additionalDistance,segment.coursePoints,segment.performance);
        
        % Taxi & deployment
        segment = FlightModel(plane,environment,missionNo,.15,0,'taxideploy',1, ...
            segment.additionalDistance,segment.coursePoints,segment.performance);
        
    else
        decelDistance = 0;
    end

    timeAtLastLapEnd = timeAtLapEnd;
    timeAtLapEnd = segment.performance.t(end);
    nLaps = nLaps + 1;
end

if timeAtLapEnd < timeLimit
    missionTime = timeAtLapEnd;
else
    missionTime = timeAtLastLapEnd;
    nLaps = nLaps - 1;
end

%% OUTPUT TO PLANETOOLS
missionResults.TOFL = TOFL;
missionResults.VTO = VTO;
missionResults.nLaps = nLaps;
missionResults.missionTime = missionTime;
missionResults.missionPoints = missionPoints;
missionResults.coursePoints = segment.coursePoints;
missionResults.performance = segment.performance;