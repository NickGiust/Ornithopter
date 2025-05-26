% Input plane, environment, mission, throttle, mission segment, additional distance,
% beginning velocity, and beginning state of charge

% Output speed vs. time, state of charge vs. time vectors and straightaway
% additional distance

function segmentResults = FlightModel(plane, environment, missionNo, throttle, CL, courseSegment, ...
    upwind, additionalDist, pointsElapsed, performanceElapsed)

% Declare global variables
global T Pel omega CTf CPf

%% Aircraft parameters
if missionNo == 2
    m = plane.m2;
    CD0 = plane.CD02;
    
    % Propeller
    D = plane.D2;

    % Ct coefficients
    A = plane.A2;
    B = plane.B2;
    C = plane.C2;

    % Cp coefficients
    E = plane.E2;
    F = plane.F2;
    G = plane.G2;
    
    % Battery
    nSeries = plane.nSeries2;
    nParallel = plane.nParallel2;
    Eb = plane.Eb2;   
    Rt = plane.Rt2;
    capacity = plane.bat2capacity;
    Rbat = plane.bat2R;
    
elseif missionNo == 3
    m = plane.m3;
    CD0 = plane.CD03;
    
    % Propeller
    D = plane.D3;

    % Ct coefficients
    A = plane.A3;
    B = plane.B3;
    C = plane.C3;

    % Cp coefficients
    E = plane.E3;
    F = plane.F3;
    G = plane.G3;
    
    % Battery
    nSeries = plane.nSeries3;
    nParallel = plane.nParallel3;
    Eb = plane.Eb3;
    Rt = plane.Rt3;
    capacity = plane.bat3capacity;
    Rbat = plane.bat3R;
end

b = plane.b;
c = plane.c;
e = plane.e;

CLmax = plane.CLmax;

rr = plane.rr;
CLgroundRoll = plane.CLgroundRoll;
hWing = plane.hWing;
pitchRate = .4;

% Motor
Kv = plane.Kv;
nMotors = plane.nMotors;

% Battery decay model constants
u = 3.5;
v = -0.0334;
w = -0.106;
z = 0.7399;
n = 1.403;
o = 2;

%% Environment parameters
windDirection = environment.windDirection;
windSpeed = environment.windSpeed;
runwayDirectionUpwind = environment.runwayDirectionUpwind;
dens = environment.dens;

%% Standard equations
%S = plane.S;
W = m*9.81;
AR = plane.AR;
S = (plane.b-plane.wFuse)*plane.c;
k = 1/(pi*e*AR);
theta = windDirection-runwayDirectionUpwind;
gfx = ((16*hWing/b)^2)/(1+((16*hWing/b)^2));
CLrate = pitchRate*2*pi/(1 + 2*pi*k);
Vstall = sqrt((2*W)/(CLmax*dens*S)); % for level flight. Multiply by load factor for accelerated flight

%% Forces on the aircraft in motion
% All forces can be split into component parts - for instance, the
% aircraft's induced drag and parasite drag in cruise.  We can rearrange Newton's 2nd
% law of motion, F = m*a, as a = F/m.  Since the aircraft, of constant
% mass, is the only object being acted upon by these forces, we can
% likewise split the aircraft's various accelerating forces (although all these 
% are decelerating) into their component parts, as follows.

% Component of aircraft deceleration due to induced drag in cruise ...
indDec = ((2*k*(W^2))/(dens*S))/m;
% ... in turns
indDecTurn = dens*S*k/(2*m);
% ... in climbs
indDecClimb = dens*S*k*(CL^2)/(2*m);
% ... due to parasite drag
paraDec = dens*S*CD0/(2*m);
% ... due to induced and parasite drag on takeoff
indParaDecTO = ( dens*S* ( CD0 + (k*gfx*(CLgroundRoll^2)) - (rr*CLgroundRoll ) ) )/(2*m);
% ... due to rolling resistance of the wheels on takeoff
wheelDecTO = rr*W/m;

% Since thrust is so complex, it's defined in the separate file
% ThrustModel.m

%% Initial values
try tInit = performanceElapsed.t(end); catch tInit = 0; end
try Vinit = performanceElapsed.V(end); catch Vinit = 0; end
try sInit = performanceElapsed.s(end); catch sInit = 1; end
try checkStall = performanceElapsed.stall(end); catch checkStall = false; end

%% Differential equation solver parameters
tmax = 50;
interval = 0.1;
tspan = tInit+interval:interval:tInit+tmax;
%% TAKEOFF
if strcmp(courseSegment, 'takeoff')
    % Takeoff speed & headwind
    segmentResults.Vstall = Vstall;
    VTO = 1.2*Vstall;
    segmentResults.VTO = VTO;
    headwind = windSpeed*cos(theta);

    % Simulate takeoff
    [t,p] = ode45(@TOEquation,0:0.1:10,[sInit,headwind]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions and trim to desired takeoff speed
    V = V(V == real(V));
    V = V(V < VTO);
    t = t(1:length(V));
    s = s(1:length(V));
    drift = 0;

    % Distance as a function of t up to takeoff speed
    x = [0; cumsum(diff(t).*(V(1:end-1)-headwind))];
   
    % Takeoff field length in FEET
    segmentResults.TOFL = x(end);
    
%% CLIMB
elseif strcmp(courseSegment, 'climb')
    % Define desired cruise altitude
    cruiseAlt = 10;
    
    % Simulate turning flight (reduced throttle and fixed CL near max)
    [t,p] = ode45(@climbEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));

    % Get drag vector
    drag = indDecClimb*(V.^2) + paraDec*(V.^2);
    
    % Get thrust vector
    thrust = nMotors*T(A,B,C,D,E,F,G,Kv,pi,Rt,V,u,v,w,z,n,nSeries,o,dens,s,throttle);
    
    % Get vertical speed (rate of climb) vector 
    Vy = ((thrust - drag).*V)/W;
    
    % Get altitude vector
    y = cumsum(Vy)*interval;
    
    % Trim to desired distance
    y = y(y < cruiseAlt);
    t = t(1:length(y));
    V = V(1:length(y));
    s = s(1:length(y));
    drift = 0;

    if sum(V < Vstall) > 0
        checkStall = true;
    end

%% STRAIGHTAWAY
elseif strcmp(courseSegment, 'halfStraightaway')
    % Standard DBF 500 ft (153 m) half-straightaway length
    flightSegmentDist = 153;
        
    % Simulate cruise flight
    [t,p] = ode45(@cruiseEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));
    
    % Adjust wind angle (theta) by 180 degrees if downwind
    if ~upwind
        theta = theta - pi;
    end

    % Calculate groundspeed
    Vg = ((windSpeed^2) + (V.^2) - (2*windSpeed*V*cos((2*pi)-theta))).^.5;

    % Calculated distance
    x = cumsum(Vg)*interval;
    
    % Trim to segment length
    x = x(x < (flightSegmentDist + additionalDist));
    t = t(1:length(x));
    V = V(1:length(x));
    s = s(1:length(x));
    drift = 0;

    if sum(V < Vstall) > 0
        checkStall = true;
    end

%% DECELERATION
elseif strcmp(courseSegment, 'deceleration')
    % Simulate decelerating cruise flight
    [t,p] = ode45(@decelEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));
    
    % Adjust wind angle (theta) by 180 degrees if downwind
    if ~upwind
        theta = theta - pi;
    end

    % Calculate groundspeed
    Vg = ((windSpeed^2) + (V.^2) - (2*windSpeed*V*cos((2*pi)-theta))).^.5;

    % Calculated distance
    x = cumsum(Vg)*interval;
    
    % Calculate landing speed (same as TO speed)
    Vland = 1.2*Vstall;
    
    % Trim to segment length - segment ends when wing approaches stall
    x = x(V > Vland);
    t = t(1:length(x));
    V = V(1:length(x));
    s = s(1:length(x));
    drift = -x(end);

    if sum(V < Vstall) > 0
        checkStall = true;
    end
    
%% STOPPING
elseif strcmp(courseSegment, 'stopping')
    % Simulate braking on the ground
    t = transpose(tInit+interval:interval:tInit+10);
    V = Vinit - 1.1*9.81*(t-tInit);
    
    x = cumsum(V)*interval;
    
    % Trim to segment length - segment ends when groundspeed is 0
    x = x(V > 0);
    t = t(1:length(x));
    V = V(1:length(x));
    s = sInit*ones(size(x));
    drift = 0;
    
%% TAXI & DEPLOY
elseif strcmp(courseSegment, 'taxideploy')    
    % Trim to segment length - segment ends when groundspeed is 0
    t = transpose(tInit+.1:.1:tInit+5.1);
    V = [ones(25,1)*2.5; zeros(26,1)];
    s = ones(length(t),1)*sInit;
    drift = 0;

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
    
    % Calculate initial lift coefficient based on initial airspeed
    CLinit = 2*W/(dens*(Vinit^2)*S);
    
    % Simulate turning flight (reduced throttle and fixed CL near max)
    [t,p] = ode45(@turnEquation,tspan,[sInit,Vinit]);
    s = p(:,1);
    V = p(:,2);

    % Eliminate nonreal portions
    V = V(V == real(V));
    t = t(1:length(V));
    s = s(1:length(V));
    
    % Get CL vector
    CL = transpose(min([CLrate*(transpose(t)-tInit)+CLinit; CL*ones(1,length(V))]));

    % Get lift vector
    L = 0.5*dens*(V.^2)*S.*CL;

    % Get load factor vector
    loadFactor = L/W;
    
    % Get r vector
    r = (V.^2)./(9.81*loadFactor);

    % Get runway-relative bearing vector
    bearing = cumsum(atan((V*interval)./r));

    % Trim to desired turn angle
    bearing = bearing(bearing < flightSegmentTurnAngle);
    loadFactor = loadFactor(1:length(bearing));
    CL = CL(1:length(bearing));
    t = t(1:length(bearing));
    V = V(1:length(bearing));
    s = s(1:length(bearing));
    VstallTime = sqrt(loadFactor).*Vstall.*ones(size(bearing));

%     if sum(V < loadFactor*Vstall) > 0
%         checkStall = true;
%     end

    % Calculate wind drift and return as additional distance for next
    % straightaway
    drift = windSpeed*(t(end)-tInit);
end

% Add wind drift during turns to results
segmentResults.additionalDistance = drift;

% Add load factor if the segment wasn't a turn
try loadFactor = loadFactor;
catch loadFactor = ones(size(t));
end

% Add level flight stall speed if segment wasn't a turn
try VstallTime = VstallTime;
catch VstallTime = Vstall * ones(size(t));
end

%% Lap Segment Markers
% Segment
segmentResults.coursePoints.segments = [pointsElapsed.segments; convertCharsToStrings(courseSegment)];
% Time
segmentResults.coursePoints.startTimes = [pointsElapsed.startTimes; tInit];

%% Time-Dependent Performance Stats
% STALL (checking if it stalled, boolean)
try segmentResults.performance.stall = [performanceElapsed.stall; checkStall];
catch segmentResults.performance.stall = checkStall; end

% STALL SPEED
try segmentResults.performance.Vstall = [performanceElapsed.Vstall; VstallTime];
catch segmentResults.performance.Vstall = VstallTime; end

% TIME VECTOR
try segmentResults.performance.t = [performanceElapsed.t; t];
catch segmentResults.performance.t = t; end

% Velocity
try segmentResults.performance.V = [performanceElapsed.V; V];
catch segmentResults.performance.V = V; end

% Battery state of charge (SoC)
try segmentResults.performance.s = [performanceElapsed.s; s];
catch segmentResults.performance.s = s; end

% Lift Coefficient
notCruiseFlight = max(strcmp(courseSegment, ... 
    ["takeoff", "taxideploy", "stopping", "climb", "deceleration"]));
levelFlight = max(strcmp(courseSegment, "halfStraightaway"));
if notCruiseFlight
    CL = ones(length(V),1)*CL;
elseif levelFlight
    CL = 2*W./(dens*(V.^2)*S);
end
try segmentResults.performance.CL = [performanceElapsed.CL; CL];
catch segmentResults.performance.CL = CL; end

% Load factor
if ~max(strcmp(courseSegment, ["turn180", "turn360"]))
    loadFactor = ones(length(V),1);
end
try segmentResults.performance.loadFactor = [performanceElapsed.loadFactor; loadFactor];
catch segmentResults.performance.loadFactor = loadFactor; end

% Induced drag power
Pinduced = .5*dens*(V.^3).*(loadFactor.^2)*k.*(CL.^2)*S;
if throttle == 0; Pinduced = zeros(size(s)); end
try segmentResults.performance.Pinduced = [performanceElapsed.Pinduced; Pinduced];
catch segmentResults.performance.Pinduced = Pinduced; end

% Parasite drag power
Pparasite = .5*dens*(V.^3)*CD0*S;
if throttle == 0; Pparasite = zeros(size(s)); end
try segmentResults.performance.Pparasite = [performanceElapsed.Pparasite; Pparasite];
catch segmentResults.performance.Pparasite = Pparasite; end

% Open-circuit voltage
Vb = nSeries*(u + v*(-log(s)).^n + w*s + z*exp(o*(s-1)));
try segmentResults.performance.Vb = [performanceElapsed.Vb; Vb];
catch segmentResults.performance.Vb = Vb; end

% Thrust
thrust = nMotors*T(A,B,C,D,E,F,G,Kv,pi,Rt,V,u,v,w,z,n,nSeries,o,dens,s,throttle);
if throttle == 0; thrust = zeros(size(s)); end
try segmentResults.performance.thrust = [performanceElapsed.thrust; thrust];
catch segmentResults.performance.thrust = thrust; end
    
% Prop speed
propspeed = omega(D,E,F,G,Kv,pi,Rt,V,u,v,w,z,n,nSeries,o,dens,s,throttle); % [rad/s]
if throttle == 0; propspeed = zeros(size(s)); end
try segmentResults.performance.propspeed = [performanceElapsed.propspeed; propspeed];
catch segmentResults.performance.propspeed = propspeed; end

% Advance ratio
J = (2*pi*V)./(propspeed*D);
try segmentResults.performance.J = [performanceElapsed.J; J];
catch segmentResults.performance.J = J; end

% Propeller coefficients
Cp = CPf*(E*(J.^2)) + (F*J) + G;
Ct = CTf*(A*(J.^2)) + (B*J) + C;
try segmentResults.performance.Cp = [performanceElapsed.Cp; Cp];
catch segmentResults.performance.Cp = Cp; end
try segmentResults.performance.Ct = [performanceElapsed.Ct; Ct];
catch segmentResults.performance.Ct = Ct; end

% Shaft power (total)
Pshaft = nMotors*Cp.*dens.*((propspeed/(2*pi)).^3)*(D^5);
if throttle == 0; Pshaft = zeros(size(s)); end
try segmentResults.performance.Pshaft = [performanceElapsed.Pshaft; Pshaft];
catch segmentResults.performance.Pshaft = Pshaft; end

% Torque
proptorque = Pshaft./(2*pi*propspeed*nMotors);
if throttle == 0; proptorque = zeros(size(s)); end
try segmentResults.performance.proptorque = [performanceElapsed.proptorque; proptorque];
catch segmentResults.performance.proptorque = proptorque; end

% Current
I = Kv*proptorque*nMotors;
try segmentResults.performance.I = [performanceElapsed.I; I];
catch segmentResults.performance.I = I; end

% Electrical power
Pelectric = I.*Vb*throttle;
try segmentResults.performance.Pelectric = [performanceElapsed.Pelectric; Pelectric];
catch segmentResults.performance.Pelectric = Pelectric; end

% Battery load (C-rate)
cRate = I/capacity;
try segmentResults.performance.cRate = [performanceElapsed.cRate; cRate];
catch segmentResults.performance.cRate = cRate; end

% Aircraft power
Paircraft = thrust.*V;
try segmentResults.performance.Paircraft = [performanceElapsed.Paircraft; Paircraft];
catch segmentResults.performance.Paircraft = Paircraft; end

% Under-load voltage
VbUnderLoad = Vb - (I*Rbat);
try segmentResults.performance.VbUnderLoad = [performanceElapsed.VbUnderLoad; VbUnderLoad];
catch segmentResults.performance.VbUnderLoad = VbUnderLoad; end

% Motor efficiency (η_motor)
etaMotor = Pshaft./Pelectric;
try segmentResults.performance.etaMotor = [performanceElapsed.etaMotor; etaMotor];
catch segmentResults.performance.etaMotor = etaMotor; end

% Propeller efficiency (η_prop)
etaProp = Paircraft./Pshaft;
try segmentResults.performance.etaProp = [performanceElapsed.etaProp; etaProp];
catch segmentResults.performance.etaProp = etaProp; end

% Overall propulsive efficiency (η)
eta = etaMotor.*etaProp;
try segmentResults.performance.eta = [performanceElapsed.eta; eta];
catch segmentResults.performance.eta = eta; end

%% Cruise system of differential equations
function dpdt = cruiseEquation(t,p)
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -indDec/((p(2))^2) - paraDec*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end

%% Turn system of differential equations
function dpdt = turnEquation(t,p)
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -(indDecTurn*min([(CLrate*(t-tInit)+CLinit)^2 CL^2]))*((p(2))^2) - paraDec*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end

%% Takeoff system of differential equations
function dpdt = TOEquation(t,p)
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);

if (plane.D2/plane.P2) >= 1.3 || (plane.D3/plane.P3) >= 1.3
    rampTime = 1;
else
    rampTime = 3;
end

if t >= 0 && t <= rampTime
    throttleTO = t * (throttle/rampTime);
elseif t > rampTime
    throttleTO = throttle;
end

% throttleTO = throttle;
Vb = nSeries*(u + v*(-log(p(1)))^n + w*p(1) + z*exp(o*(p(1)-1)));

% Coefficients to solve for propeller speed with quadratic formula after
% setting motor torque equal to propeller torque
Ca = CPf*(G*dens*(D^5))/(8*(pi^3));
Cb = CPf*((F*dens*p(2)*(D^4))/(4*(pi^2))) - (1/((Kv^2)*Rt));
Cc = CPf*((E*dens*(p(2)^2)*(D^3))/(2*pi)) + ((throttleTO*Vb)/(Kv*Rt));

omegatest = (-Cb - ( ( (Cb^2) - (4*Ca*Cc) ).^(1/2) )) / (2*Ca);

J = (2*pi*p(2))/(omegatest*D);

% % T = CTf * dens * (p(2)^2) * (D^2) * (A + (B/J) + (C/(J^2)))
% % T = CTf * dens * (p(2)^2) * (D^2) * (A*J^2 + B*J + C)/(J^2)
% % T = CTf * dens * (p(2)^2) * (D^2) * (A*J^2 + B*J + C)/(((2*pi*p(2))/(omega*D))^2)
% % T = CTf * dens * (p(2)^2) * (D^2) * (A*J^2 + B*J + C)*(omega^2*D^2)/(4 * pi^2 * p(2)^2)
Ttest = CTf * dens * (omegatest^2) * (D^4) * (A*J^2 + B*J + C) / (4 * pi^2);

Cp = CPf*((E*(J^2)) + (F*J) + G);

N = omegatest/(2*pi);
Ps = Cp*(N^3)*dens*(D^5);

Peltest = Vb*throttleTO*Ps*Kv/omegatest;

if p(2) > 0
    dVdt = -indParaDecTO*((p(2))^2) - wheelDecTO + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttleTO);
else
    dVdt = -indParaDecTO*((p(2))^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttleTO);
end
dpdt = [dsdt; dVdt];

end

%% Climb system of differential equations
function dpdt = climbEquation(t,p)
dsdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dVdt = -indDecClimb*((p(2))^2) - paraDec*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dsdt; dVdt];
end


%% Deceleration system of differential equations
function dpdt = decelEquation(t,p)
dsdt = 0;
dVdt = -indDec/((p(2))^2) - paraDec*(p(2)^2);
dpdt = [dsdt; dVdt];
end

end