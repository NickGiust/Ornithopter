clear;
clc;

% wind direction
windDirection = 117*pi/180;

% wind speed
windSpeed = 1;

% runway direction
runwayDirectionUpwind = 105*pi/180;

theta = windDirection-runwayDirectionUpwind;

% Input plane, mission, throttle, course segment, beginning velocity, and beginning state of
% charge
% Output speed vs. time, state of charge vs. time vectors and straightaway
% additional distance

global A B C D E F G Kv Rt u v w z n nSeries o dens throttle M Mturn N nMotors m Eb

% Airframe
m = .795;

b = 1.88;
c = .28;
e = .864;

CD0 = .0313;
CLmax = 1.11;

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
Kv = 40.5;
Rt = .18;
nMotors = 1;

% Battery
% 3S lithium polymer pack
throttle = 1;
nSeries = 3;
Eb = 2.2*3.7*3600;
u = 3.5;
v = -.0334;
w = -.106;
z = .74;
n = 1.4;
o = 2;
                                                                                        
% Airfield
dens = 1.225;

S = b*c;
W = m*9.81;
Vstall = sqrt((2*W)/(CLmax*dens*S));
VTO = 1.2*Vstall;
AR = (b^2)/S;
k = 1/(pi*e*AR);

M = ((2*k*(W^2))/(dens*S))/m;
CLturn = .9*CLmax;
Mturn = dens*S*(CLturn^2)/(2*m);
N = ((dens*S*CD0)/2)/m;

%% Calculation parameters
tmax = 50;
interval = 0.1;
tspan = 0:interval:tmax;

%% 1000 ft (305 m) STRAIGHTAWAY
% Simulate cruise flight
[t,p] = ode45(@cruiseEquation,tspan,[1,8]);
s = p(:,1);
V = p(:,2);

% Eliminate nonreal portions
V = V(V == real(V));
t = t(1:length(V));
s = s(1:length(V));

% Calculate groundspeed
Vg = (V./sin(pi-theta)) .* (theta - ((windSpeed*sin(pi-theta))./V));

% Trim to distance
x = cumsum(Vg)*interval;
x = x(x < 305);
t = t(1:length(x));
V = V(1:length(x));
s = s(1:length(x));


%% 500 ft (158 m) HALF-STRAIGHTAWAY
% Simulate cruise flight

% Trim to distance

%% 360-DEGREE MAX LIFT TURN
% Simulate turning flight (reduced throttle and fixed CL near max)
[t,p] = ode45(@turnEquation,tspan,[1,4]);
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

% Trim to 360
bearing = bearing(bearing < 2*pi);
t = t(1:length(bearing));
V = V(1:length(bearing));
s = s(1:length(bearing));

% Calculate groundspeed magnitude


% Calculate groundspeed components

% Calculate x and y in ground reference frame

% Calculate x and y in air reference frame

% Take magnitude of difference and return as additional distance for next
% straightaway

plot(t,V)
hold on;
plot(t,s*10);
hold off;

%% 180-DEGREE MAX LIFT TURN
% Simulate cruise flight (reduced throttle and fixed CL near max)

% Get lift vector

% Get load factor vector

% Get r vector

% Get runway-relative bearing vector

% Trim to 360

% Calculate groundspeed magnitude

% Calculate groundspeed components

% Calculate x and y in ground reference frame

% Calculate x and y in air reference frame

% Take magnitude of difference and return as additional distance for next
% straightaway

function dpdt = cruiseEquation(t,p)
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle M N nMotors m Eb

%% Thrust with respect to state of charge
T = @(A,B,C,D,E,F,G,Kv,PI,Rt,V,a,b,c,d,n,nSeries,o,rho,s,throttle) ...
    A.*D.^2.*V.^2.*rho+(C.*1.0./D.^6.*1.0./G.^2.*PI.^4.*...
    (-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0)...
    .^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)...
    -(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))...
    ./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)...
    ./4.0).^2.*4.0)./rho-(B.*1.0./D.^2.*PI.^2.*V.*(-sqrt((1.0./...
    Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./...
    PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*...
    (a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0)...
    +1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).*2.0)./G;

%% Electrical power with respect to state of charge
Pel = @(D,E,F,G,Kv,PI,Rt,V,a,b,c,d,n,nSeries,o,rho,s,throttle)...
    (1.0./D.^5.*1.0./G.^2.*Kv.*PI.^3.*nSeries.*throttle.*(G-(D.^4.*F...
    .*G.*1.0./PI.^3.*V.*rho.*pi)./(sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*...
    1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*...
    V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+...
    b.*(-log(s)).^n))./(Kv.*Rt)))./2.0).*-2.0+(1.0./Kv.^2.*2.0)./Rt+...
    (D.^4.*F.*1.0./PI.^2.*V.*rho)./2.0)+(D.^8.*E.*G.^2.*1.0./PI.^6.*V.^2.*...
    rho.^2.*pi.^2.*1.0./(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*...
    V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./...
    (PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*...
    (-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^...
    2.*V.*rho)./4.0).^2)./4.0).*(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./...
    PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*...
    rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*...
    (-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^...
    2.*V.*rho)./4.0).^2.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n).*2.0)./rho;

%% Cruise system of differential equations
dVdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dsdt = -M/((p(2))^2) - N*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dVdt; dsdt];
end


function dpdt = turnEquation(t,p)
global A B C D E F G Kv Rt u v w z n nSeries o dens throttle Mturn N nMotors m Eb

%% Thrust with respect to state of charge
T = @(A,B,C,D,E,F,G,Kv,PI,Rt,V,a,b,c,d,n,nSeries,o,rho,s,throttle) ...
    A.*D.^2.*V.^2.*rho+(C.*1.0./D.^6.*1.0./G.^2.*PI.^4.*...
    (-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0)...
    .^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)...
    -(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))...
    ./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)...
    ./4.0).^2.*4.0)./rho-(B.*1.0./D.^2.*PI.^2.*V.*(-sqrt((1.0./...
    Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./...
    PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*...
    (a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0)...
    +1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).*2.0)./G;

%% Electrical power with respect to state of charge
Pel = @(D,E,F,G,Kv,PI,Rt,V,a,b,c,d,n,nSeries,o,rho,s,throttle)...
    (1.0./D.^5.*1.0./G.^2.*Kv.*PI.^3.*nSeries.*throttle.*(G-(D.^4.*F...
    .*G.*1.0./PI.^3.*V.*rho.*pi)./(sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*...
    1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*...
    V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+...
    b.*(-log(s)).^n))./(Kv.*Rt)))./2.0).*-2.0+(1.0./Kv.^2.*2.0)./Rt+...
    (D.^4.*F.*1.0./PI.^2.*V.*rho)./2.0)+(D.^8.*E.*G.^2.*1.0./PI.^6.*V.^2.*...
    rho.^2.*pi.^2.*1.0./(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*...
    V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./...
    (PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*...
    (-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^...
    2.*V.*rho)./4.0).^2)./4.0).*(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./...
    PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*...
    rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*...
    (-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^...
    2.*V.*rho)./4.0).^2.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n).*2.0)./rho;

%% Turn system of differential equations
dVdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dsdt = -Mturn*((p(2))^2) - N*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dVdt; dsdt];
end