clc; clear; close all;

%% Inputs
currentDirectory = split(pwd, 'PlaneTools22');
planeToolsDirectory = [currentDirectory{1} '\PlaneTools22'];
airplane = 'Prototype1TF2';
year = '2022S';

%% Read aircraft from text file
plane = ReadPlane(airplane, year, planeToolsDirectory);

%% Trade Study

% chordArray = (plane.c*0.5):0.01:(plane.c*2); % chord TS
% payloadArray = (plane.mPayload2*1):0.1:(plane.mPayload2*2); % payload TS
% scoreArray = zeros(length(chordArray),length(payloadArray));
% 
% figure
% for i = 1:length(chordArray)
%     for j = 1:length(payloadArray)
%         plane.c = chordArray(i);
%         plane.mPayload2 = payloadArray(j);
%         results = TS_SimulateCompetition(plane, planeToolsDirectory);
%         scoreArray(i,j) = results.score;
%         results.TOFL
%     end
% end
% 
% [X,Y] = meshgrid(chordArray,payloadArray);
% surf(X,Y,scoreArray')
% 
% xlabel("chord [m]")
% ylabel("payload [kg]")
% zlabel("score")

results = TS_SimulateCompetition(plane, planeToolsDirectory);