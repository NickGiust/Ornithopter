clc; clear; close all;
global planeToolsDirectory;

%% SIMULATION CORE PARAMETER
currentDirectory = split(pwd, 'PlaneTools22');
planeToolsDirectory = [currentDirectory{1} '\PlaneTools22'];
year = '2022S';
airplane = 'Prototype1TF3';
location = 'Wichita';

%% Add directories to path and initialize global propulsion functions
addpath([planeToolsDirectory '\PhysicsModels\']);
addpath([planeToolsDirectory '\ComponentLibrary']);
addpath([planeToolsDirectory '\Plotting']);
addpath([planeToolsDirectory '\Environments']);
global T CTf CPf
CTf = 1;
CPf = 1;
if ~isa(T,'function_handle') 
    ThrustModel; 
end

%% Read aircraft from text file
plane = ReadPlane(airplane, year, planeToolsDirectory);

%% Read environment from text file 
environment = ReadEnvironment(location, planeToolsDirectory);

%% MISSION 2
% Build up weight and drag
plane.e = dragBuildup(plane, 2).e;
plane.m2 = massBuildup(plane, 2).m;
%plane.m2 = 3.09;
plane.CD02 = dragBuildup(plane, 2).CD0;

% Set mission parameters
mission.missionNo = 2;
mission.lapLimit = 40;
mission.timeLimit = 600;
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.crzThrottle = 0.9;
mission.turnThrottle = 0.8;
mission.climbCL = .6;
mission.turnCL = 2.1;

% chord trade study
chordArray = (plane.c*0.5):0.01:(plane.c*1.5);
chordArray = plane.c;
lapsArray = [];
% Simulate mission
for i = 1:length(chordArray)
    plane.c = chordArray(i);
    missionResults2 = SimulateMission(plane, environment, mission);
    lapsArray = [lapsArray missionResults2.halfLoadedLaps];
end

figure
plot(chordArray,lapsArray)
xlabel("chord [m]")
ylabel("mission laps")

PlotPerformance(mission, plane, missionResults2);
%CheckLimits(mission, plane, missionResults2);

% Assumptions
% Score mission


%% Competition scoring

