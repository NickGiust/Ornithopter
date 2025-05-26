clc;
clear;

Dinches = 14;
Pinches = 7;
prefix = 'apce_';

propName = [prefix num2str(Dinches) 'x' num2str(Pinches)];
D = Dinches*.0254;
dens = 1.225;

[~,cmdout] = system(['qprop Propellers/' propName '.txt Motors/hacker_a60_6xs_v4_28_pole.txt 0,80/30 200,15000/30 0']);
rawData = splitlines(convertCharsToStrings(cmdout));
rawData = rawData(18:end,:);

propdata = [];
for r = 1:size(rawData,1)
    numData = split(rawData(r));
    numData = numData(numData ~= "");
    numData = str2double(numData);
    try propdata = [propdata; transpose(numData)];
    end
end

V = propdata(:,1);
omegaRevPerMin = propdata(:,2);
T = propdata(:,4);
P = propdata(:,6);

n = omegaRevPerMin/60;

J = V./(n*D);
CT = T./(dens*(n.^2)*(D^4));
CP = P./(dens*(n.^3)*(D^5));

combined = [J CT CP];
combined = combined((combined(:,2) > 0 & combined(:,3) > 0),:);
J = combined(:,1);
CT = combined(:,2);
CP = combined(:,3);

scatter(J,CT)
ylim([0 .12])

CTfit = fit(J,CT,'poly2');
ABC = coeffvalues(CTfit);

CPfit = fit(J,CP,'poly2');
EFG = coeffvalues(CPfit);

propeller.D = D;
propeller.A = ABC(1);
propeller.B = ABC(2);
propeller.C = ABC(3);
propeller.E = EFG(1);
propeller.F = EFG(2);
propeller.G = EFG(3);

eval([propName ' =  propeller; clear propeller'])
eval(['clearvars -except ' propName]);
%load('propeller.mat');
%global planeToolsDirectory;
%save([planeToolsDirectory '\ComponentLibrary\propeller.mat']);