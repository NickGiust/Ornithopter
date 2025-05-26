clc; clear; close all;
global planeToolsDirectory;

%% SIMULATION CORE PARAMETER
currentDirectory = split(pwd, 'PlaneTools22');
planeToolsDirectory = [currentDirectory{1} '\PlaneTools22'];
year = '2022';
airplane = 'Prototype1TF2';
location = 'Wichita';

%% Add directories to path and initialize global propulsion functions
addpath([planeToolsDirectory '\PhysicsModels\']);
addpath([planeToolsDirectory '\ComponentLibrary']);
addpath([planeToolsDirectory '\Plotting']);
addpath([planeToolsDirectory '\Environments']);
global T CTf CPf
CTf = 1;
CPf = 1;
if ~isa(T,'function_handle'); ThrustModel; end

%% Read aircraft from text file
plane = ReadPlane(airplane, year, planeToolsDirectory);

%% Read environment from text file 
environment = ReadEnvironment(location, planeToolsDirectory);

%% MISSION 2
% Build up weight and drag
plane.e = dragBuildup(plane, 3).e;
plane.m2 = massBuildup(plane, 2).m;
plane.CD02 = dragBuildup(plane, 2).CD0;

% Set mission parameters
mission.missionNo = 2;
mission.lapLimit = 3;
mission.timeLimit = 10000;
mission.TOthrottle = 1;
mission.climbThrottle = 1;
mission.crzThrottle = 1;
mission.turnThrottle = .9;
mission.climbCL = .6;
mission.turnCL = 2.1;

% Simulate mission
missionResults2 = SimulateMission(plane, environment, mission);

PlotPerformance(mission, plane, missionResults2);
%CheckLimits(mission, plane, missionResults2);

% Assumptions
% Score mission

%% MISSION 3
% Build up weight and drag
plane.m3 = massBuildup(plane, 3).m;
plane.CD03 = dragBuildup(plane, 3).CD0;

% Set mission parameters
mission.missionNo = 3;
mission.lapLimit = 100;
mission.timeLimit = 600;
mission.TOthrottle = .9;
mission.climbThrottle = 1;
mission.crzThrottle = 1;
mission.turnThrottle = 1;
mission.climbCL = .6;
mission.turnCL = 2.1;

% Simulate mission
missionResults3 = SimulateMission(plane, environment, mission);

%PlotPerformance(mission, plane, missionResults3);
%CheckLimits(mission, plane, missionResults3);

% Assumptions
% Score mission

%% Competition scoring
