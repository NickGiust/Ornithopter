syms V PI omega D

J = (2*PI*V)/(omega*D);

syms X Y Z P

etaP = X*((J/P)^2) + Y*(J/P) + Z;

omega = solve(etaP, omega);
omega = omega(1);

syms OMEGA throttle Vb Kv Rt Ps

Ps = (throttle*Vb)/(Kv*Rt);


Ca = -1/((Kv^2)*Rt);
Cb = (throttle*Vb)/(Kv*Rt);
Cc = -Ps;

OMEGA = (-Cb+(((Cb^2)-(4*Ca*Cc))^(1/2)))/(2*Ca);