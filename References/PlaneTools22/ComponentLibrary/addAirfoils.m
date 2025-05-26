clc;
clear;

global planeToolsDirectory;

e560.CLmax2D = 2.1;
e560.thicknessRatio = .1;
e560.wingThicknessLocation = .2;

ba527ls.CLmax2D = 1.51;
ba527ls.thicknessRatio = .1;
ba527ls.wingThicknessLocation = .2;

e423.CLmax2D = 2.3;
e423.thicknessRatio = .15;
e423.wingThicknessLocation = .2;

e422.CLmax2D = 2.1;
e422.thicknessRatio = .15;
e422.wingThicknessLocation = .24;

save([planeToolsDirectory '\ComponentLibrary\airfoil.mat']);