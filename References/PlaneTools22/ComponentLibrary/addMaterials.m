clc;
clear;


balsa.density = 160; %kg/m3
balsa.roughness = .0005;

plywood.density = 680; %kg/m3

carbonFiber.density = 5.5 / 29.494; %oz/yd2 --> kg/yd2

solite.density = 0.7 / 29.494; %oz/yd2 --> kg/yd2


% Archive
steel.density = 8000;
steel.roughness = 0;

carbon.density = 1.5;
carbon.roughness = .001;

foamfiberglass.density = 2.05;
foamfiberglass.roughness = .004;

kevlar_plywood.density = 70;
kevlar_plywood.roughness = .0009;

% global planeToolsDirectory;
currentDirectory = split(pwd, 'PlaneTools22');
planeToolsDirectory = [currentDirectory{1} '\PlaneTools22'];
save([planeToolsDirectory '\ComponentLibrary\material.mat']);
