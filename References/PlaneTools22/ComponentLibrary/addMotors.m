clc;
% clear;

% for Kv, use registered Kv/(2pi/60)

Scorpion_HKIII_5025_520KV_F3S.Kv = 54.454;
Scorpion_HKIII_5025_520KV_F3S.R = .009;
Scorpion_HKIII_5025_520KV_F3S.m = .595;
Scorpion_HKIII_5025_520KV_F3S.maxPower = 7696;
Scorpion_HKIII_5025_520KV_F3S.length = 0.0649;
Scorpion_HKIII_5025_520KV_F3S.diam = 0.0615;

Hacker_A50_8S_V3.Kv = 89;
Hacker_A50_8S_V3.R = .03;
Hacker_A50_8S_V3.m = .26;
Hacker_A50_8S_V3.maxPower = 2800;

Hacker_A60_5S_V4_28pole.Kv = 31;
Hacker_A60_5S_V4_28pole.R = .015;
Hacker_A60_5S_V4_28pole.m = .6;
Hacker_A60_5S_V4_28pole.maxPower = 2600;

Hacker_A60_6XS_V4_28pole.Kv = 38.7;
Hacker_A60_6XS_V4_28pole.R = .015;
Hacker_A60_6XS_V4_28pole.m = .6;
Hacker_A60_6XS_V4_28pole.maxPower = 2300;

Hacker_A50_10L_Turnado_V3.Kv = 55.5;
Hacker_A50_10L_Turnado_V3.R = .015;
Hacker_A50_10L_Turnado_V3.m = .445;
Hacker_A50_10L_Turnado_V3.maxPower = 6000;

Scorpion_MIII_4015_310kv.Kv = 310*(2*pi)/60;
Scorpion_MIII_4015_310kv.R = .109;
Scorpion_MIII_4015_310kv.m = .208;
Scorpion_MIII_4015_310kv.maxPower = 962;

Hacker_B50_15XL.Kv = 1052*(2*pi)/60;
Hacker_B50_15XL.R = .279;
Hacker_B50_15XL.m = .34;
Hacker_B50_15XL.maxPower = 1700;

Turnigy_Aerodrive_SK3_5045_660KV.Kv = 660*(2*pi)/60;
Turnigy_Aerodrive_SK3_5045_660KV.R = .014;
Turnigy_Aerodrive_SK3_5045_660KV.m = .280;
Turnigy_Aerodrive_SK3_5045_660KV.maxPower = 1410;

Hacker_A50_12L_V3.Kv = 450*(2*pi)/60;
Hacker_A50_12L_V3.R = .019;
Hacker_A50_12L_V3.m = .445;
Hacker_A50_12L_V3.maxPower = 4500;
Hacker_A50_12L_V3.length = 0.061;
Hacker_A50_12L_V3.diam = 0.0515;

Neu_Motors_1915_1Y.Kv = 590*(2*pi)/60;
Neu_Motors_1915_1Y.R = 0.015;
Neu_Motors_1915_1Y.m = 0.396;
Neu_Motors_1915_1Y.maxPower = 3570;

Scorpion_SII_4020_540KV.Kv = 540*(2*pi)/60;
Scorpion_SII_4020_540KV.R = 0.02;
Scorpion_SII_4020_540KV.m = 0.288;
Scorpion_SII_4020_540KV.maxPower = 1850;
Scorpion_SII_4020_540KV.length = 0.04615; % m
Scorpion_SII_4020_540KV.diam = 0.0489; % m

TMotor_AT4125_540Kv.Kv = 540*(2*pi)/60;
TMotor_AT4125_540Kv.R = 0.014;
TMotor_AT4125_540Kv.m = 0.423; % 0.355
TMotor_AT4125_540Kv.maxPower = 2000;
TMotor_AT4125_540Kv.length = 0.050;
TMotor_AT4125_540Kv.diam = 0.050;

Scorpion_A_4220_540.Kv = 540*(2*pi)/60;
Scorpion_A_4220_540.R = 0.0182;
Scorpion_A_4220_540.m = 0.288;
Scorpion_A_4220_540.maxPower = 2553;
Scorpion_A_4220_540.length = 0.0578;
Scorpion_A_4220_540.diam = 0.0507;

%global planeToolsDirectory;
save([planeToolsDirectory '\ComponentLibrary\motor.mat']);

