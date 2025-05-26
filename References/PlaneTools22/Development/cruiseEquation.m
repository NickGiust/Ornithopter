function dpdt = cruiseEquation(t,p)

% Airframe
m = .795;

b = 1.88;
c = .28;
e = .864;

CD0 = .0313;
CLmax = 1.5;

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
% 3S lithium polymer pack
throttle = .3;
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
N = ((dens*S*CD0)/2)/m;

%% Thrust with respect to state of charge
T = @(A,B,C,D,E,F,G,Kv,PI,Rt,V,a,b,c,d,n,nSeries,o,rho,s,throttle)A.*D.^2.*V.^2.*rho+(C.*1.0./D.^6.*1.0./G.^2.*PI.^4.*(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2.*4.0)./rho-(B.*1.0./D.^2.*PI.^2.*V.*(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).*2.0)./G;

%% Electrical power with respect to state of charge
Pel = @(D,E,F,G,Kv,PI,Rt,V,a,b,c,d,n,nSeries,o,rho,s,throttle)(1.0./D.^5.*1.0./G.^2.*Kv.*PI.^3.*nSeries.*throttle.*(G-(D.^4.*F.*G.*1.0./PI.^3.*V.*rho.*pi)./(sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0).*-2.0+(1.0./Kv.^2.*2.0)./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./2.0)+(D.^8.*E.*G.^2.*1.0./PI.^6.*V.^2.*rho.^2.*pi.^2.*1.0./(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2)./4.0).*(-sqrt((1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2-(D.^5.*G.*1.0./PI.^3.*rho.*((D.^3.*E.*V.^2.*rho)./(PI.*2.0)-(nSeries.*throttle.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n))./(Kv.*Rt)))./2.0)+1.0./Kv.^2./Rt+(D.^4.*F.*1.0./PI.^2.*V.*rho)./4.0).^2.*(a+c.*s+d.*exp(o.*(s-1.0))+b.*(-log(s)).^n).*2.0)./rho;

%% Cruise system of differential equations
dVdt = (1/Eb) * -Pel(D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dsdt = -M/((p(2))^2) - N*(p(2)^2) + ...
        (nMotors/m)*T(A,B,C,D,E,F,G,Kv,pi,Rt,p(2),u,v,w,z,n,nSeries,o,dens,p(1),throttle);
dpdt = [dVdt; dsdt];
end