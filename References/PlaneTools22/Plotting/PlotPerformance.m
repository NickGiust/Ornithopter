function PlotPerformance(mission, plane, missionResults)
performance = missionResults.performance;

%% TIME-BASED PLOTS
% Airspeed
figure;
subplot(2,1,1)
plot(performance.t,performance.V,'Color','black','LineWidth',2)
ylabel('Airspeed [m/s]');
ylim([0 max(performance.V)*1.3]);
title([' Mission ' num2str(mission.missionNo) ' Performance']);
plotLapDivisions;

% State of charge
subplot(2,1,2)
plot(performance.t,performance.s*100,'Color','black','LineWidth',2)
ylabel('State of Charge [%]');
xlabel('Time [s]');
ylim([min(performance.s)*100 min(performance.s)*100+(max(performance.s) ...
    -min(performance.s))*120]);
plotLapDivisions;

% Prop speed
figure;
subplot(3,1,1)
plot(performance.t,performance.propspeed*60/(2*pi),'Color','black','LineWidth',1.5);
ylabel('Prop. speed [rev/min]');
ylim([0 max(performance.propspeed)*1.3*60/(2*pi)]);
plotLapDivisions;

% Current/battery load
subplot(3,1,2)
plot(performance.t,performance.I,'Color','black','LineWidth',1.5);
ylabel('Current [A]');
ylim([0 max(performance.I)*1.3]);
yyaxis right;
plot(performance.t,performance.cRate,'Color','black','LineWidth',1.5);
ylabel('C-rate [1/h]');
ylim([0 max(performance.cRate)*1.3]);
plotLapDivisions;

% Voltage
subplot(3,1,3)
hold on;
plot(performance.t,performance.VbUnderLoad,'Color','red','LineWidth',1.5);
plot(performance.t,performance.Vb,'Color','green','LineWidth',1.5);
hold off;
ylabel('Battery voltage [V]');
xlabel('Time [s]');
plotLapDivisions;

% Lift coefficient
figure;
subplot(2,1,1)
plot(performance.t,performance.CL,'Color','black','LineWidth',1.5);
ylabel('Lift coefficient');
xlabel('Time [s]');
ylim([0 2.5]);
plotLapDivisions;

% Load factor
subplot(2,1,2);
plot(performance.t,performance.loadFactor,'Color','black','LineWidth',1.5);
ylabel('Load factor, n');
xlabel('Time [s]');
plotLapDivisions;

% Stacked line chart of power to battery, motor, prop, induced, parasite drag
figure;
dragPower = performance.Pinduced + performance.Pparasite;
inducedShare = performance.Pinduced./dragPower;
parasiteShare = performance.Pparasite./dragPower;

area(performance.t, [performance.Paircraft.*inducedShare ...
    performance.Paircraft.*parasiteShare ...
    performance.Pshaft - performance.Paircraft ...
    performance.I.^2 * (plane.Rt2 - (plane.bat2R/plane.nParallel2)) ...
    performance.I.^2 * (plane.bat2R/plane.nParallel2)]);
ylabel('Power [W]');
xlabel('Time [s]');
legend("Induced drag", "Parasite drag", "Propeller", "Motor/ESC/Wires", "Battery");

end