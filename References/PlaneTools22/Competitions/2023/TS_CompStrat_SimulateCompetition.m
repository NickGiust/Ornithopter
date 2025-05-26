function [missionResults2, missionResults3, missionResultsGM, plane] = TS_CompStrat_SimulateCompetition(plane, environment, planeToolsDirectory)

% Why this script? No mass convergence loop!

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

% Build up weight and drag
[M2massResults, plane] = massBuildup(plane, 2);

plane.m2 = M2massResults.m;
plane.trueGMloadingMultiplier = checkGMTiming(plane);
missionResultsGM.score = 650/plane.m2;

%% MISSION 2
% Build up weight and drag
dragResults2 = dragBuildup(plane, 2);
plane.CD02 = dragResults2.CD0;
plane.e = dragBuildup(plane, 2).e;
plane.dragResults2 = dragResults2;
plane.CD02 = dragBuildup(plane, 2).CD0;

% Set mission parameters
mission.missionNo = 2;
mission.lapLimit = 40;
mission.timeLimit = 600; % 10 minutes
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.climbCL = 1.0;
mission.turnCL = 0.85*plane.CLmax;

TOFLLimit = 60*0.3048;
nLapsArray = [];
missionResults2Array = [];
crzThrottleList = 0.3:0.05:0.6;

for i = 1:length(crzThrottleList)
    mission.crzThrottle = crzThrottleList(i);
    mission.turnThrottle = mission.crzThrottle;
    TOthrottleConvergence = true;
    missionResults2 = SimulateMission(plane, environment, mission, TOFLLimit, TOthrottleConvergence);
    nLapsArray = [nLapsArray missionResults2.nLaps];
    missionResults2Array = [missionResults2Array missionResults2];
end

[~, maxLapsIndex] = max(nLapsArray);

missionResults2 = missionResults2Array(maxLapsIndex);
missionResults2.optimizedThrottleSetting = crzThrottleList(maxLapsIndex);

if missionResults2.minTOthrottle < 1 && missionResults2.TOFLbool == "PASSED"
    TOthrottleConvergence = false;
    missionResults2_fullTOthrottle = SimulateMission(plane, environment, mission, TOFLLimit,TOthrottleConvergence);
    plane.m2TOFL_fullThrottle = missionResults2_fullTOthrottle.TOFL;
    TOsegmentIndex = find(missionResults2_fullTOthrottle.performance.V > missionResults2_fullTOthrottle.VTO, 1);
    plane.m2TOFL_avgCurrent = mean(missionResults2_fullTOthrottle.performance.I(1:TOsegmentIndex-1));
    plane.m2TOFL_avgCurrentTime = missionResults2_fullTOthrottle.coursePoints.startTimes(2);
else
    plane.m2TOFL_fullThrottle = 0;
    plane.m2TOFL_avgCurrent = 0;
    plane.m2TOFL_avgCurrentTime = 0;
end

% % Plot throttle vs. nLaps
% figure
% for i = 1:length(missionResults2Array)
%     hold on
%     scatter(crzThrottleList(i), missionResults2Array(i).nLaps,'blue','filled')
%     xlabel('Throttle Setting')
%     ylabel('Number of Laps')
%     title('M2 Throttle Optimization')
%     ax = gca;
%     ax.FontSize = 18;
% end



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
mission.TOthrottle = 0.95;
mission.climbThrottle = 0.95;
mission.crzThrottle = 0.95;
mission.turnThrottle = 0.95;
mission.climbCL = 0.6;
mission.turnCL = 0.8*plane.CLmax;

if missionResults2.TOFLbool == "FAILED"
    missionResults3.score = 0;
    missionResults3.missionNo = 3;
    missionResults3.TOFLbool = "NOT ATTEMPTED";
    missionResults3.TOFL = 0;
    missionResults3.nLaps = 0;
    missionResults3.missionTime = 0;
    missionResults3.minTOthrottle = 1.1;
    plane.m3TOFL_fullThrottle = 0;
    plane.m3TOFL_avgCurrent = 0;
    plane.m3TOFL_avgCurrentTime = 0;
else
    % Simulate mission
    TOthrottleConvergence = true;
    missionResults3 = SimulateMission(plane, environment, mission, TOFLLimit,TOthrottleConvergence);
end

if missionResults3.minTOthrottle < 1 && missionResults3.TOFLbool == "PASSED"
    TOthrottleConvergence = false;
    missionResults3_fullTOthrottle = SimulateMission(plane, environment, mission, TOFLLimit,TOthrottleConvergence);
    plane.m3TOFL_fullThrottle = missionResults3_fullTOthrottle.TOFL;
    TOsegmentIndex = find(missionResults3_fullTOthrottle.performance.V > missionResults3_fullTOthrottle.VTO, 1);
    plane.m3TOFL_avgCurrent = mean(missionResults3_fullTOthrottle.performance.I(1:TOsegmentIndex-1));
    plane.m3TOFL_avgCurrentTime = missionResults3_fullTOthrottle.coursePoints.startTimes(2);
end

% PlotPerformance(mission, plane, missionResults3);
% CheckLimits(mission, plane, missionResults3);





end

