clc;
clear;
global planeToolsDirectory;

Turnigy_LiPo_4500.capacity = 4.5;
Turnigy_LiPo_4500.R = .03;
Turnigy_LiPo_4500.Imax = 300;
Turnigy_LiPo_4500.mPerCell = .13;

% 6s or 4s
HRB_LiPo_3300.capacity = 3.3;
HRB_LiPo_3300.R = .033;
HRB_LiPo_3300.Imax = 330;
HRB_LiPo_3300.mPerCell = .089;
HRB_LiPo_3300.nSeries = 8;

% For 10S
Thunder_LiPo_2700.capacity = 2.7;
Thunder_LiPo_2700.R = .036;
Thunder_LiPo_2700.Imax = 231;
Thunder_LiPo_2700.mPerCell = .687/8;

% For 8S
Thunder_LiPo_3300.capacity = 3.3;
Thunder_LiPo_3300.R = .036;
Thunder_LiPo_3300.Imax = 231;
Thunder_LiPo_3300.mPerCell = .687/8;

% For 6S
Thunder_LiPo_4400.capacity = 4.4;
Thunder_LiPo_4400.R = .036;
Thunder_LiPo_4400.Imax = 308;
Thunder_LiPo_4400.mPerCell = .691/6;

% For 4S
Thunder_LiPo_6600.capacity = 6.6;
Thunder_LiPo_6600.R = .036;
Thunder_LiPo_6600.Imax = 462;
Thunder_LiPo_6600.mPerCell = .676/4;

% % TP1800-6SR70
% Thunder_LiPo_1800.capacity = 1.8;
% % Thunder_LiPo_1800.R = ;
% Thunder_LiPo_1800.Imax = 252;
% Thunder_LiPo_1800.mass = 0.291;
% Thunder_LiPo_1800.nSeries = 6;
% 
% % HRB-Power 1800 100C
% HRB_LiPo_1800.capacity = 1.8;
% % HRB_LiPo_1800.R = ;
% % HRB_LiPo_1800.Imax = ;
% HRB_LiPo_1800.mass = 0.270;
% HRB_LiPo_1800.nSeries = 6;

%% New batteries from https://www.thunderpowerrc.com/collections/elite-series-55c

Thunder_1300_6S.capacity = 1.3;
Thunder_1300_6s.R = 0.036;
Thunder_1300_6s.Imax = 71.5;
Thunder_1300_6s.mass = 0.200;
Thunder_1300_6s.nSeries = 6;

Thunder_2250_6s.capacity = 2.25;
Thunder_2250_6s.R = 0.036;
Thunder_2250_6s.Imax = 123.8;
Thunder_2250_6s.mass = 0.345;
Thunder_2250_6s.nSeries = 6;

Thunder_2700_6s.capacity = 2.7;
Thunder_2700_6s.R = 0.036;
Thunder_2700_6s.Imax = 148.5;
Thunder_2700_6s.mass = 0.409;
Thunder_2700_6s.nSeries = 6;

Thunder_3300_6s.capacity = 3.3;
Thunder_3300_6s.R = 0.036;
Thunder_3300_6s.Imax = 181.5;
Thunder_3300_6s.mass = 0.510;
Thunder_3300_6s.nSeries = 6;

Thunder_3300_7s.capacity = 3.3;
Thunder_3300_7s.R = 0.036;
Thunder_3300_7s.Imax = 181.5;
Thunder_3300_7s.mass = 0.590;
Thunder_3300_7s.nSeries = 7;

Thunder_3300_8s.capacity = 3.3;
Thunder_3300_8s.R = 0.036;
Thunder_3300_8s.Imax = 181.5;
Thunder_3300_8s.mass = 0.674;
Thunder_3300_8s.nSeries = 8;

Thunder_3300_9s.capacity = 3.3;
Thunder_3300_9s.R = 0.036;
Thunder_3300_9s.Imax = 181.5;
Thunder_3300_9s.mass = 0.758;
Thunder_3300_9s.nSeries = 9;

save([planeToolsDirectory '\ComponentLibrary\battery.mat']);