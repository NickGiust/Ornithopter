clear;
clc;

%% CRUISE CALCULATIONS
% Calculation parameters
tmax = 50;
interval = 0.1;

% V as a function of t up to cruise speed
tspan = 0:interval:tmax;

[t,p] = ode45(@cruiseEquation,tspan,[1,5]);
s = p(:,1);
V = p(:,2);
V = V(V == real(V));
t = t(1:length(V));
s = s(1:length(V));

plot(t,V)
hold on;
plot(t,s*10);
hold off;

%% Statistics AT 50% OF FLIGHT
%{
Vcruise = 
Vb = 

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

%}