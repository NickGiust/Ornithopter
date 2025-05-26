% Input plane, environment, mission, throttle, mission segment, additional distance,
% beginning velocity, and beginning state of charge

% Output speed vs. time, state of charge vs. time vectors and straightaway
% additional distance

function segmentPerformance = FlightModel(plane, environment, missionNo, throttle, CL, courseSegment, Vinit)

% Declare global variables
global T omega

%% Aircraft parameters
if missionNo == 3 || missonNo == 2
    m = plane.m3;
    CD0 = plane.CD03;
    
    % Propeller
    D = plane.D3;
    
    % Pitch
    P = plane.P;
    
    CPv = [-1.0239    5.3990    0.0383   -8.9173   -0.4666    0.0696    4.8682    0.7180   -0.2624];
    CTv = [-0.7917    4.3418   -0.0606   -6.5257   -1.1431    0.3439    3.1725    1.5924   -0.5897];

    % Ct coefficients
    A = CTv(6) + (CTv(9) * P);
    B = CTv(3) + (CTv(5) * P) + (CTv(8) * (P^2));
    C = CTv(1) + (CTv(2) * P) + (CTv(4) * (P^2)) + (CTv(7) * (P^3));

    % Cp coefficients
    E = CPv(6) + (CPv(9) * P);
    F = CPv(3) + (CPv(5) * P) + (CPv(8) * (P^2));
    G = CPv(1) + (CPv(2) * P) + (CPv(4) * (P^2)) + (CPv(7) * (P^3));
    
    % Battery
    Vb = plane.Vb;
    Rt = plane.Rt3;
    batMaxCurrent = plane.bat3maxCurrent*plane.nParallel3;
end

b = plane.b;
c = plane.c;
e = plane.e;

CLmax = plane.CLmax;

rr = plane.rr;
CLgroundRoll = plane.CLgroundRoll;
hWing = plane.hWing;
nStruct = plane.nStruct;

% Motor
Kv = plane.Kv;
nMotors = plane.nMotors;
maxPowerAllMotors = plane.motorMaxPower*plane.nMotors;

% ESC
ESCMaxCurrent = plane.ESCMaxCurrent;

% General η for propeller stall
propEtaMin = .1;

% Advance ratio cutoff for quadratic CP and CT model
Jmin = .25;

%% Environment parameters
dens = environment.dens;

%% Standard equations
S = b*c;
W = m*9.81;
AR = (b^2)/S;
k = 1/(pi*e*AR);
M = ((2*k*(W^2))/(dens*S))/m;
Mturn = dens*S*k*(CL^2)/(2*m);
N = dens*S*CD0/(2*m);
gfx = ((16*hWing/b)^2)/(1+((16*hWing/b)^2));
MTO = (dens*S*(CD0 + (k*gfx*(CLgroundRoll^2)) - (rr*CLgroundRoll)))/(2*m);
NTO = rr*W/m;

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
    [t,p] = ode45(@TOEquation,tspan,0);
    V = p;

    % Eliminate nonreal portions and trim to desired takeoff speed
    V = V(V == real(V));
    V = V(V < VTO);
    t = t(1:length(V));
    
    % Distance as a function of t up to takeoff speed
    x = cumsum(diff(t).*V(1:end-1));

    % Takeoff field length in FEET
    segmentPerformance.TOFL = x(end);

%% STRAIGHTAWAY
elseif strcmp(courseSegment, 'halfStraightaway')
    % Standard DBF 500 ft (153 m) half-straightaway length
    flightSegmentDist = 153;
        
    % Simulate cruise flight
    [t,p] = ode45(@cruiseEquation,tspan,20);
    V = p;
    
    % Cruise wing stall warning - ensures that CL < CLmax
    CL = 2*W./(dens*(V.^2)*S);
    if max(CL) > CLmax
        %cruiseCL = max(CL)
        %CLmax
        error('Wing stall warning: CL exceeds wing maximum 3D CL.')
    end
    
%% MAXIMUM-LIFT TURN
else
    % Standard DBF 180-degree turn at each end of the course
    flightSegmentTurnAngle = pi;
    
    % Simulate turning flight (reduced throttle and fixed CL near max)
    [t,p] = ode45(@turnEquation,tspan,Vinit);
    V = p;

    % Get lift vector
    L = 0.5*dens*(V.^2)*S*CL;

    % Get load factor vector
    loadFactor = L/W;
    
    % Exceeding turn structural limit warning - ensures that n < nStruct
    if max(loadFactor) > nStruct
        %nTurn = max(loadFactor)
        %nStruct
        error('Structure warning: aircraft structural load factor limit exceeded during turn.');
    end
    
    % Get r vector
    r = (V.^2)./(9.81*loadFactor);

    % Get runway-relative bearing vector
    bearing = cumsum(atan((V*interval)./r));

    % Trim to desired turn angle
    bearing = bearing(bearing < flightSegmentTurnAngle);
    t = t(1:length(bearing));
    V = V(1:length(bearing));
    
end

segmentPerformance.t = t(end);
segmentPerformance.V = V(end);

%% Propulsion Warnings
thrust = nMotors*T(A,B,C,D,E,F,G,Kv,pi,Rt,V,Vb,dens,throttle);
propspeed = omega(D,E,F,G,Kv,pi,Rt,V,Vb,dens,throttle);

J = (2*pi*V)./(propspeed*D);
Cp = (E*(J.^2)) + (F*J) + G;
Pshaft = Cp.*dens.*((propspeed/(2*pi)).^3)*(D^5);
proptorque = Pshaft./propspeed;
I = Kv*proptorque;
Pelectric = I.*Vb*throttle;
Paircraft = thrust.*V/nMotors;
etaMotor = Pshaft./Pelectric;
etaProp = Paircraft./Pshaft;

segmentPerformance.Pel = Pelectric(end);

%% Cruise system of differential equations
function dpdt = cruiseEquation(t,p)
dVdt = -M/(p^2) - N*(p^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p,Vb,dens,throttle);
dpdt = [dVdt];
end


%% Turn system of differential equations
function dpdt = turnEquation(t,p)
dVdt = -Mturn*(p^2) - N*(p^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p,Vb,dens,throttle);
dpdt = [dVdt];
end

%% Takeoff system of differential equations
function dpdt = TOEquation(t,p)
dVdt = -MTO*(p^2) - NTO + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p,Vb,dens,throttle);
dpdt = [dVdt];
end

end