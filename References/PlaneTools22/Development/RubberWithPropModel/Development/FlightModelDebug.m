clear;
clc;

% Input plane, environment, mission, throttle, mission segment, additional distance,
% beginning velocity, and beginning state of charge

% Output speed vs. time, state of charge vs. time vectors and straightaway
% additional distance

% Declare global variables
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle M Mturn N ...
    nMotors m Eb MTO NTO Mclimb T

%% Aircraft parameters (***ALL VALUES IN SI BASE UNITS***)
m = 9.6;

b = 1.22;
c = .41;
e = .8;

CD0 = .04;
CLmax = 1.2;

rr = .07;
CLgroundRoll = .4;
hWing = .2;

CLclimb = .6;

% Propeller
% APC 17x12E
D = 17*.0254;
propVals = [-0.1413   -0.0092    0.1047   -0.1669    0.1150    0.0293];

% Ct coefficients
A = propVals(1);
B = propVals(2);
C = propVals(3);

% Cp coefficients
E = propVals(4);
F = propVals(5);
G = propVals(6);

% Motor
% Hacker A40-14L-V4 355 Kv
Kv = 37.2;
Rt = .09;
nMotors = 2;

% Battery
% 6S2P Turnigy LiPo 4500mAh
nSeries = 6;
Eb = 9*3.7*3600;
u = 3.5;
v = -.0334;
w = -.106;
z = .74;
n = 1.4;
o = 2;

%% Environment parameters
windDirection = 117*pi/180;
windSpeed = 1;
runwayDirectionUpwind = 105*pi/180;
dens = 1.225;

%% DIRECT FUNCTION INPUTS*****
mission = 2;

throttle = .9;

courseSegment = 'takeoff';
%courseSegment = 'climb';
%courseSegment = 'straightaway';
%courseSegment = 'turn180';
%courseSegment = 'halfStraightaway';
%courseSegment = 'turn360';

additionalDist = 5;

Vinit = 15;

sInit = 1;

%% Standard equations
S = b*c;
W = m*9.81;
Vstall = sqrt((2*W)/(CLmax*dens*S));
VTO = 1.2*Vstall;
AR = (b^2)/S;
k = 1/(pi*e*AR);
theta = windDirection-runwayDirectionUpwind;
M = ((2*k*(W^2))/(dens*S))/m;
CLturn = .9*CLmax;
Mturn = dens*S*(CLturn^2)/(2*m);
N = ((dens*S*CD0)/2)/m;
gfx = ((16*hWing/b)^2)/(1+((16*hWing/b)^2));
MTO = (dens*S*(CD0 + (k*gfx*(CLgroundRoll^2)) - (rr*CLgroundRoll)))/(2*m);
NTO = (rr*W)/m;
Mclimb = dens*S*(CLclimb^2)/(2*m);

%% Differential equation solver parameters
tmax = 50;
interval = 0.1;
tspan = 0:interval:tmax;

%% TAKEOFF
if strcmp(courseSegment, 'takeoff')
    % Takeoff speed
    Vstall = sqrt((2*W)/(CLmax*dens*S));
    VTO = 1.2*Vstall;

    % Simulate takeoff
    [t,p] = ode45(@TOEquation,tspan,[1,0]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions and trim to desired takeoff speed
    V = V(V == real(V));
    V = V(V < VTO);
    t = t(1:length(V));
    s = s(1:length(V));
    
    subplot(2,1,1)
    plot(t,V)
    ylabel('Airspeed [m/s]');

    subplot(2,1,2)
    plot(t,s*100)
    ylabel('State of Charge [%]');
    xlabel('Time [s]');
    
%% CLIMB
elseif strcmp(courseSegment, 'climb')
    % Define desired cruise altitude
    cruiseAlt = 15;
    
    % Simulate turning flight (reduced throttle and fixed CL near max)
    [t,p] = ode45(@climbEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));

    % Get drag vector
    drag = Mclimb*(V.^2) + N*(V.^2);
    
    % Get thrust vector
    thrust = nMotors*T(A,B,C,D,E,F,G,Kv,pi,Rt,V,u,v,w,z,n,nSeries,o,dens,s,throttle);
    
    % Get RoC vector 
    RoC = ((thrust - drag).*V)/W;
    
    % Get altitude vector
    y = cumsum(RoC)*interval;
    
    % Trim to desired distance
    y = y(y < cruiseAlt);
    t = t(1:length(y));
    V = V(1:length(y));
    s = s(1:length(y));
    
    subplot(2,1,1)
    plot(t,V)
    ylabel('Airspeed [m/s]');

    subplot(2,1,2)
    plot(t,s*100)
    ylabel('State of Charge [%]');
    xlabel('Time [s]');

%% STRAIGHTAWAY
elseif strcmp(courseSegment, 'straightaway') || strcmp(courseSegment, 'halfStraightaway')
    if strcmp(courseSegment, 'straightaway')
        % Standard DBF 1000 ft (305 m) straightaway length
        flightSegmentDist = 305;
    else
        % Standard DBF 500 ft (153 m) half-straightaway length
        flightSegmentDist = 153;
    end
        
    % Simulate cruise flight
    [t,p] = ode45(@cruiseEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));

    % Calculate groundspeed
    Vg = (V./sin(pi-theta)) .* (theta - ((windSpeed*sin(pi-theta))./V));

    % Trim to desired distance
    x = cumsum(Vg)*interval;
    x = x(x < flightSegmentDist);
    t = t(1:length(x));
    V = V(1:length(x));
    s = s(1:length(x));
    
    subplot(2,1,1)
    plot(t,V)
    ylabel('Airspeed [m/s]');

    subplot(2,1,2)
    plot(t,s*100)
    ylabel('State of Charge [%]');
    xlabel('Time [s]');

%% MAXIMUM-LIFT TURN
else
    if strcmp(courseSegment, 'turn180')
        % Standard DBF 180-degree turn at each end of the course
        flightSegmentTurnAngle = pi;
    else
        % Standard DBF 360-degree turn at the midpoint of downwind
        % straightaway
        flightSegmentTurnAngle = 2*pi;
    end
    % Simulate turning flight (reduced throttle and fixed CL near max)
    [t,p] = ode45(@turnEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));

    % Get lift vector
    L = 0.5*dens*(V.^2)*S*CLturn;

    % Get load factor vector
    loadFactor = L/W;

    % Get r vector
    r = (V.^2)./(9.81*loadFactor);

    % Get runway-relative bearing vector
    bearing = cumsum(atan((V*interval)./r));

    % Trim to desired turn angle
    bearing = bearing(bearing < flightSegmentTurnAngle);
    t = t(1:length(bearing));
    V = V(1:length(bearing));
    s = s(1:length(bearing));

    % Calculate wind drift and return as additional distance for next
    % straightaway
    drift = windSpeed*t(end);

    subplot(2,1,1)
    plot(t,V)
    ylabel('Airspeed [m/s]');

    subplot(2,1,2)
    plot(t,s*100)
    ylabel('State of Charge [%]');
    xlabel('Time [s]');
end


%% Cruise system of differential equations
function dpdt = cruiseEquation(t,p)
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle M N nMotors m Eb T Pel
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -M/((p(2))^2) - N*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end

%% Turn system of differential equations
function dpdt = turnEquation(t,p)
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle Mturn N nMotors m Eb T Pel
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -Mturn*((p(2))^2) - N*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end

%% Takeoff system of differential equations
function dpdt = TOEquation(t,p)
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle MTO NTO nMotors m Eb T Pel
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -MTO*((p(2))^2) - NTO + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end

%% Climb system of differential equations
function dpdt = climbEquation(t,p)
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle Mclimb N nMotors m Eb T Pel
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -Mclimb*((p(2))^2) - N*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end