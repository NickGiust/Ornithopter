clear;
clc;

%% PARAMETERS - ****ALL VALUES IN SI BASE UNITS****

% Airframe
% AF0010 new airframe

%***CHANGE***
config = 'Gliding Power Consumption - Config. A - Battery in Nacelle - ';

%***CHANGE***
m = .795;

b = 1.88;
c = .28;
e = .864;

%***CHANGE***
CD0 = .0313;

CLmax = 1.5;

%***CHANGE***
DAshares = [0.7491    0.0524    0.0198    0.0216    0.1363    0.0207];
DAlabels = {'Wing', 'Nacelle', 'Motor', 'Gaps', 'Winglets', 'Propeller'};

%***CHANGE***
mShares = [0.3635    0.2566    0.0252    0.0503    0.0377    0.0503    0.0629    0.0252    0.0566    0.0126    0.0277  0.0314];
mLabels = {'Airframe (W)', 'Battery (W)', 'Elevons (W)', 'ESC (W)', 'Flight Controller (W)', ...
    'Motor (W)', 'Nacelle (W)', 'Pitot (W)', 'Propeller (W)', 'Receiver (W)', 'Servos (W)', ...
    'Stick/Motor Mount (W)'};

%***CHANGE***
explode = [0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 1 0 0];

% Propeller
% Aeronaut CAM 14x9
D = 14*.0254;
propVals = [-0.1422   -0.0507    0.1342   -0.1905    0.1294    0.0298];

% Ct coefficients
A = propVals(1);
B = propVals(2);
C = propVals(3);

% Cp coefficients
E = propVals(4);
F = propVals(5);
G = propVals(6);

% Motor
% Hacker A10-9L with 4.4:1 gearbox
Kv = 178;
Rt = .18;
nMotors = 1;

% Battery
% 3S lithium polymer pack, with 8% loss
Vcell = 3.7;
Nseries = 3;
Vloss = .08;
% Minimum Power
%throttle = .169;
% Minimum Drag
throttle = .177
                                                                                        
% Airfield
dens = 1.225;

% Calculation parameters
tmax = 50;
interval = 0.1;

S = b*c;
W = m*9.81;
Vstall = sqrt((2*W)/(CLmax*dens*S));
VTO = 1.2*Vstall;
AR = (b^2)/S;
k = 1/(pi*e*AR);
Vb = Vcell*Nseries*(1-Vloss);

%% Dynamic thrust equation
T = @(A,B,C,D,E,F,G,Kv,PI,Rt,V,Vb,rho,throttle) A.*D.^2.*V.^2.*rho+...
    (C.*1.0./D.^6.*1.0./G.^2.*PI.^4.*...
    (-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2+...
    (D.^5.*G.*1.0./PI.^3.*rho.*((Vb.*throttle)./...
    (Kv.*Rt)-(D.^3.*E.*V.^2.*rho)./(PI.*2.0)))./2.0)+1.0./...
    Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2.*4.0)./rho-...
    (B.*1.0./D.^2.*PI.^2.*V.*(-sqrt((1.0./Kv.^2./Rt+...
    (D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2+(D.^5.*G.*1.0./PI.^3.*...
    rho.*((Vb.*throttle)./(Kv.*Rt)-(D.^3.*E.*V.^2.*rho)./(PI.*2.0)))./2.0)+...
    1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).*2.0)./G;

%% CRUISE CALCULATIONS
% V as a function of t up to cruise speed
tspan = 0:interval:tmax;
M = ((2*k*(W^2))/(dens*S))/m;
N = ((dens*S*CD0)/2)/m;
[t,V] = ode45((@(t,V) (T(A,B,C,D,E,F,G,Kv,pi,Rt,V,Vb,dens,throttle)*nMotors/m)...
    -(M/(V^2))-(N*(V^2))), tspan, VTO);

Vcruise = V(end)

%% Statistics
CL = 2*W/(dens*(Vcruise^2)*S)
induced = ((2*k*(W^2))/(dens*S))/(Vcruise^2)
parasite = ((dens*S*CD0)/2)*(Vcruise^2)
thrust = T(A,B,C,D,E,F,G,Kv,pi,Rt,Vcruise,Vb,dens,throttle)*nMotors
liftToDrag = W/(induced + parasite)

omega = @(D,E,F,G,Kv,PI,Rt,V,Vb,rho,throttle)(1.0./D.^5.*PI.^3.*...
    (-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2+...
    (D.^5.*G.*1.0./PI.^3.*rho.*((Vb.*throttle)./(Kv.*Rt)-...
    (D.^3.*E.*V.^2.*rho)./(PI.*2.0)))./2.0)+1.0./Kv.^2./Rt+...
    (D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).*-4.0)./(G.*rho);
propspeed = omega(D,E,F,G,Kv,pi,Rt,Vcruise,Vb,dens,throttle)

J = (2*pi*Vcruise)/(propspeed*D);
Cp = (E*(J^2)) + (F*J) + G;
Pshaft = Cp*dens*((propspeed/(2*pi))^3)*(D^5);
proptorque = Pshaft/propspeed;
I = Kv*proptorque
Pelectric = I*Vb*throttle;
Paircraft = thrust*Vcruise/nMotors;
etaMotor = Pshaft/Pelectric;
etaProp = Paircraft/Pshaft;

Pelectric
etaMotor
Pshaft
etaProp
Paircraft

Vyglide = Paircraft/W

plot(t,V)
xlabel('Time [s]');
ylabel('Airspeed [m/s]');


%% Power consumption breakdown
Pshares = [(induced/thrust)*mShares (parasite/thrust)*DAshares];
Plabels = [mLabels DAlabels];

figure(2);
pie2(Pshares, explode, Plabels);
title([config ...
    'Power Req. = ' num2str(Paircraft, 3) ' W, Sink Rate = ' num2str(Vyglide, 3) ...
    ' m/s']);
legend(Plabels);
