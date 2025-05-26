close all;
wind = 0;
dens = 1.19;
m = 12.8;
W = m*9.81;
CLmax = 2.1;
S = .4*2.44;
VTOfactor = 1;
propType = ["14x7", "15x8", "16x10"];
pitch = [14 15 16];
F = [15.5 20 23]*4.45;
%F = [3.8 20 23]*4.45;

a = F/m;
VTO = VTOfactor*(W/(.5*dens*CLmax*S))^.5;
TOFLpredicted = (VTO-wind)^2./(2*a)
TOFLactual = [110 45 28]*.3048

scatter(F, TOFLpredicted)
hold on
scatter(F, TOFLactual)
legend("kinematic prediction", "actual")
xlabel("thrust [N]")
ylabel("TOFL [m]")

flight12tas = readmatrix('C:\ADT\PlaneTools22\Development\airspeed22622\Flight 12 TAS.csv');
tflight2acc = flight12tas(flight12tas(:,1) > 1534 & flight12tas(:,1) < 1542 & flight12tas(:,2) > 2, 1);
tflight2acc = tflight2acc - tflight2acc(1);
Vflight2acc = flight12tas(flight12tas(:,1) > 1534 & flight12tas(:,1) < 1542 & flight12tas(:,2) > 2, 2);

figure;
scatter(tflight2acc,Vflight2acc)
hold on;

[accfit, gof] = fit(tflight2acc,Vflight2acc,'poly1');
plot(accfit);
fitrsquare = gof.rsquare
acc = coeffvalues(accfit);
Fcalc = m*acc(1)
xlabel("t [s]")
ylabel("Airspeed, V [m/s]")
title("T/O acceleration - 14x7 propeller")



flight3tas = readmatrix('C:\ADT\PlaneTools22\Development\airspeed22622\Flight 3 TAS.csv');

tflight3acc = flight3tas(flight3tas(:,1) > 508.5 & flight3tas(:,1) < 511 & flight3tas(:,2) > 2, 1);
tflight3acc = tflight3acc - tflight3acc(1);
Vflight3acc = flight3tas(flight3tas(:,1) > 508.5 & flight3tas(:,1) < 511 & flight3tas(:,2) > 2, 2);

figure;
scatter(tflight3acc,Vflight3acc)
hold on;

[accfit, gof] = fit(tflight3acc,Vflight3acc,'poly1');
plot(accfit);
fitrsquare = gof.rsquare
acc3 = coeffvalues(accfit);
Fcalc3 = m*acc3(1)
xlabel("t [s]")
ylabel("Airspeed, V [m/s]")
title("T/O acceleration - 16x10 propeller")
