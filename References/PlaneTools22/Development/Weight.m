clear;
%% Weight Buildup Script for AF0010 Aircraft
% Jackson Markow
% 6/8/2021

%% PARAMETERS - ****ALL VALUES IN SI BASE UNITS****

% Airframe
mAirframe = 289;

% Motor
mMotor = .04;

% Propeller
mProp = .045;

% Flight controller
mFC = .03;

% Battery
mBat = .204;

% Speed controller
mESC = .04;

% Reciever
mReciever = .01;

% Stick/motor mount
mStick = .025;

% Nacelle
%mNacelle = .05;
mNacelle = 0

% Elevons
mElevons = .02;

% Pitot tube
mPitot = .02;

% Servos/horns
mServos = .022;

%% Total
m = [mAirframe mBat mElevons mESC mFC mMotor mNacelle mPitot mProp mReciever mServos mStick];
mTotal = sum(m);
mShares = m/mTotal

mLabels = {'Airframe (W)', 'Battery (W)', 'Elevons (W)', 'ESC (W)', 'Flight Controller (W)', ...
    'Motor (W)', 'Nacelle (W)', 'Pitot (W)', 'Propeller (W)', 'Receiver (W)', 'Servos (W)', ...
    'Stick/Motor Mount (W)'};