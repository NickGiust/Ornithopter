dens = 1.225;
V = 0;
s = .5;
nSeries = 6;
capacity = 4.4;
nMotors = 3;
throttle = 1;
propType = 'ancf_15x10';
motorType = 'Hacker_A50_10L_Turnado_V3';
Resc = .003;
Rwire = .005;
Rbat = .008;

planeToolsDirectory = 'C:\ADT\PlaneTools22';
addpath([planeToolsDirectory '\ComponentLibrary']);

load('propeller.mat', propType);
propeller = eval(propType);

D = propeller.D;
A = propeller.A;
B = propeller.B;
C = propeller.C;
E = propeller.E;
F = propeller.F;
G = propeller.G;

load('motor.mat', motorType);
motor = eval(motorType);

Kv = motor.Kv;
Rmotor = motor.R;

u = 3.5;
v = -.0334;
w = -.106;
z = .74;
n = 1.4;
o = 2;

Rt = Resc + Rwire + Rmotor + Rbat;

Vb = nSeries*(u + v*(-log(s)).^n + w*s + z*exp(o*(s-1)));
thrust = nMotors*T(A,B,C,D,E,F,G,Kv,pi,Rt,V,u,v,w,z,n,nSeries,o,dens,s,throttle);
propspeed = omega(D,E,F,G,Kv,pi,Rt,V,u,v,w,z,n,nSeries,o,dens,s,throttle);
J = (2*pi*V)./(propspeed*D);
Cp = (E*(J.^2)) + (F*J) + G;
Ct = (A*(J.^2)) + (B*J) + C;
Pshaft = nMotors*Cp.*dens.*((propspeed/(2*pi)).^3)*(D^5);
proptorque = Pshaft./(propspeed*nMotors);
I = Kv*proptorque*nMotors;
Pelectric = I.*Vb*throttle;
Paircraft = thrust.*V/nMotors;
etaMotor = Pshaft./Pelectric;
etaProp = Paircraft./Pshaft;
eta = etaMotor.*etaProp;

batteryLoad = Pelectric/(3.7*nSeries*capacity);
VbUnderLoad = Vb - (I*Rbat);
motorPowerLoss = Pelectric-Pshaft;

excelData = [Vb; VbUnderLoad; thrust; I; propspeed*60/(2*pi); propspeed/(2*pi); ...
    propspeed; Pelectric; batteryLoad; Cp; Ct; etaMotor; motorPowerLoss];