clear;
clc;

global omega T

%% Propeller angular velocity/prop. speed [omega]
% Propeller constants
syms D A B C E F G
% Air density
syms rho
% Pi is set as a variable to prevent MATLAB from combining it into terrible
% constants
syms PI
% Motor & resistance characteristics
syms Kv Rt 
% Throttle setting
syms throttle
% Battery constants
syms a b c d n o nSeries;
% Airspeed
syms V
% Battery voltage
syms Vb

% Coefficients to solve for propeller speed with quadratic formula after
% setting motor torque equal to propeller torque
Ca = (G*rho*(D^5))/(8*(PI^3));
Cb = ((F*rho*V*(D^4))/(4*(PI^2))) + (1/((Kv^2)*Rt));
Cc = ((E*rho*(V^2)*(D^3))/(2*PI)) - ((throttle*Vb)/(Kv*Rt));

omega = (-Cb+(((Cb^2)-(4*Ca*Cc))^(1/2)))/(2*Ca);

%% Propeller dynamic thrust [T]
T = (A*rho*(V^2)*(D^2)) + ((B*rho*V*(D^3))/(2*PI))*omega ...
    + ((C*rho*(D^4))/(4*(PI^2)))*(omega^2);

%% Convert to MATLAB functions
T = matlabFunction(T);
omega = matlabFunction(omega);