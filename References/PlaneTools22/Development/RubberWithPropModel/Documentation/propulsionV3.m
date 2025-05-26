clear
syms a b c

syms A B C
syms E F G
syms rho V D PI Kv Rt throttle Vb
syms OMEGA T

a = (G*rho*(D^5))/(8*(PI^3));
b = ((F*rho*V*(D^4))/(4*(PI^2))) + (1/((Kv^2)*Rt));
c = ((E*rho*(V^2)*(D^3))/(2*PI)) - ((throttle*Vb)/(Kv*Rt));

OMEGA = (-b+(((b^2)-(4*a*c))^(1/2)))/(2*a);
T = (A*rho*(V^2)*(D^2)) + ((B*rho*V*(D^3))/(2*PI))*OMEGA ...
    + ((C*rho*(D^4))/(4*(PI^2)))*(OMEGA^2);

omega = matlabFunction(OMEGA);
T = matlabFunction(T);

a = -.1085;
b = -.07802;
c = .1519;
D = .229;
e = -.1369;
f = .04544;
g = .07185;
Kv = 89;
Rt = .0072;
%rho = 1.13;
rho = 1.225;
%throttle = 1;
throttle = .71;
%V = 38.4;
V = 27.2;
Vb = 20.4;

speed = omega(D,e,f,g,Kv,pi,Rt,V,Vb,rho,throttle)
thrust = T(a,b,c,D,e,f,g,Kv,pi,Rt,V,Vb,rho,throttle)

syms x;
fplot(T(a,b,c,D,e,f,g,Kv,pi,Rt,x,Vb,rho,throttle))
xlim([0,60]);
xlabel('Airspeed [m/s]');
ylabel('Dynamic Thrust [N]')
title('APC 9x6E, \rho = 1.13 kg/m^3, 6S LiPo (8% loss), Hacker A50 8S V3');




J = (2*pi*V)/(speed*D)
Cp = (e*(J^2)) + (f*J) + g
Pshaft = Cp*rho*((speed/(2*pi))^3)*(D^5)
proptorque = Pshaft/speed
I = Kv*proptorque
Pelectric = I*Vb*throttle
Phalfaircraft = thrust*V
