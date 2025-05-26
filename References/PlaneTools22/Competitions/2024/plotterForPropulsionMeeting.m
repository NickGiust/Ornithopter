close all;

% This is just a script to plot some of the results of the best plane

%% M2 data

% time
m2time = bestPlane{1,1}.m2results.performance.t; % [s]

% velocity
m2V = bestPlane{1,1}.m2results.performance.V; % [m/s]
m2Vstall = bestPlane{1,1}.m2results.performance.Vstall; % [m/s]

% battery power
m2bat = bestPlane{1,1}.m2results.performance.s; % [%]

%% M2 plotting

% V vs t
figure(1)
hold on
plot(m2time, m2V, 'LineWidth', 1.5, 'Color', '#990000')
plot(m2time, m2Vstall, 'LineWidth', 1.5, 'Color', '#FFC72C', 'LineStyle', '--')
% adding in segment markers
for i = 3:length(bestPlane{1,1}.m2results.coursePoints.startTimes)
        if strcmpi(bestPlane{1,1}.m2results.coursePoints.segments(i), 'halfStraightaway')
            xline(bestPlane{1,1}.m2results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Straightaway')
        elseif strcmpi(bestPlane{1,1}.m2results.coursePoints.segments(i), 'turn180')
            xline(bestPlane{1,1}.m2results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Turn')
        elseif strcmpi(bestPlane{1,1}.m2results.coursePoints.segments(i), 'turn360')
            xline(bestPlane{1,1}.m2results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', '360')
        elseif strcmpi(bestPlane{1,1}.m2results.coursePoints.segments(i), 'takeoff')
            xline(bestPlane{1,1}.m2results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Takeoff')
        elseif strcmpi(bestPlane{1,1}.m2results.coursePoints.segments(i), 'climb')
            xline(bestPlane{1,1}.m2results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Climb')
        end
end
% formatting
title('M2 Velocity as a Function of Time')
legend('Airspeed', 'Stall Speed')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
xlim([min(m2time), max(m2time)])

% s vs t
figure(2)
hold on
plot(m2time, 100*m2bat, 'LineWidth', 1.5, 'Color', '#990000')
title('M2 Battery % as a Function of Time')
xlabel('Time [s]')
ylabel('Percent Charge [%]')
xlim([min(m2time), max(m2time)])
ylim([0, 100])

%% M3 data

% time
m3time = bestPlane{1,1}.m3results.performance.t; % [s]

% velocity
m3V = bestPlane{1,1}.m3results.performance.V; % [m/s]
m3Vstall = bestPlane{1,1}.m3results.performance.Vstall; % [m/s]

% battery power
m3bat = bestPlane{1,1}.m3results.performance.s; % [%]

%% M3 Plotting

% V vs t
figure(3)
hold on
plot(m3time, m3V, 'LineWidth', 1.5, 'Color', '#990000')
plot(m3time, m3Vstall, 'LineWidth', 1.5, 'Color', '#FFC72C', 'LineStyle', '--')
% adding in segment markers
for i = 3:length(bestPlane{1,1}.m3results.coursePoints.startTimes)
        if strcmpi(bestPlane{1,1}.m3results.coursePoints.segments(i), 'halfStraightaway')
            xline(bestPlane{1,1}.m3results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Straightaway')
        elseif strcmpi(bestPlane{1,1}.m3results.coursePoints.segments(i), 'turn180')
            xline(bestPlane{1,1}.m3results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Turn')
        elseif strcmpi(bestPlane{1,1}.m3results.coursePoints.segments(i), 'turn360')
            xline(bestPlane{1,1}.m3results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', '360')
        elseif strcmpi(bestPlane{1,1}.m3results.coursePoints.segments(i), 'takeoff')
            xline(bestPlane{1,1}.m3results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Takeoff')
        elseif strcmpi(bestPlane{1,1}.m3results.coursePoints.segments(i), 'climb')
            xline(bestPlane{1,1}.m3results.coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Climb')
        end
end
% formatting
title('M3 Velocity as a Function of Time')
legend('Airspeed', 'Stall Speed')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
xlim([min(m3time), max(m3time)])

% s vs t
figure(4)
hold on
plot(m3time, 100*m3bat, 'LineWidth', 1.5, 'Color', '#990000')
title('M3 Battery % as a Function of Time')
xlabel('Time [s]')
ylabel('Percent Charge [%]')
xlim([min(m3time), max(m3time)])
ylim([0, 100])



