clc; clear all; close all;
global planeToolsDirectory slash;

%% Inputs
currentDirectory = split(pwd, 'PlaneTools22');
currentDirectoryCharacterVector = char(currentDirectory(1,:));
backslashLocations = strfind(currentDirectoryCharacterVector, '\');
if isempty(backslashLocations)
    slash = '/'; % forward slash for Apple Mac directories
else
    slash = '\'; % back slash for Windows directories
end

planeToolsDirectory = [currentDirectory{1} 'PlaneTools22'];

%% Add directories to path and initialize global propulsion functions
addpath([planeToolsDirectory slash 'PhysicsModels' slash]);
addpath([planeToolsDirectory slash 'PhysicsModels' slash 'MassBuildup']);
addpath([planeToolsDirectory slash 'ComponentLibrary']);
addpath([planeToolsDirectory slash 'Plotting']);
addpath([planeToolsDirectory slash 'Environments']);
% addpath([planeToolsDirectory slash 'ComponentLibrary' slash 'material.mat']);

%% Read aircraft from text file
plane = definePlane();

%% Define Environment
environment = defineEnvironment("Lucerne");

%% Trade Study Setup

% Define parameter ranges and values
% % overall
% span_range = 0.50 : 0.25 : 1.50;  % [m] max is 1.52 m (5ft)
% AR_range = 3 : 1 : 8;             % [-]
% % M2
% cabinetMass_range = 0.25 : 0.25 : 2; % [kg]
% % M3
% nPassengers_range = 5: 5 : 30; % [-] only whole numbers

% overall
span_range = 1.68;  % [m] max is 1.52 m (5ft)
AR_range = 4.5;             % [-]
% M2
cabinetMass_range = 1.3; % [kg]
% M3
nPassengers_range = 3; % [-] only whole numbers

% Generate grids of parameter values
[span, AR, cabinetMass, nPassengers] = ndgrid(span_range, AR_range, cabinetMass_range, nPassengers_range);

% Generate grids of indices
[span_index, AR_index, cabinetMass_index, nPassengers_index] = ndgrid(1:numel(span_range), 1:numel(AR_range), 1:numel(cabinetMass_range), 1:numel(nPassengers_range));

% Combine parameter grids and indices
combined = cat(2, reshape(span, [], 1), reshape(AR, [], 1), reshape(cabinetMass, [], 1), reshape(nPassengers, [], 1), ...
                   reshape(span_index, [], 1), reshape(AR_index, [], 1), reshape(cabinetMass_index, [], 1), reshape(nPassengers_index, [], 1));
lengthCombined = size(combined, 1); % calculating the length once (rather than every time it's looped)

% Create a matrix to store the calculated results
planes = cell(1, lengthCombined);
TOFL = zeros(1, lengthCombined);

% Set up a timer
tic
usedTimer = 0;
timeNow = datetime('now');

%% Trade study

% Iterate over combinations
for i = 1:lengthCombined
    % Extract parameter values and indices for the current combination
    current_span = combined(i, 1);
    current_AR = combined(i, 2);
    current_cabinetMass = combined(i, 3);
    current_nPassengers = combined(i, 4);
%     current_span_index = combined(i, 4);        % not used
%     current_AR_index2 = combined(i, 5);         % not used
%     current_radarMass_index3 = combined(i, 6);  % not used

    % Print the current plane specs so user knows it's running
    fprintf('Plane:              %.0f of %.0f\nCabinet mass [kg]:  %.2f\nNum Passengers [-]: %.2f\nSpan [m]:           %.1f\nAR:                 %.1f\n\n',...
        i, lengthCombined, current_cabinetMass, current_nPassengers, current_span, current_AR);
    
    % Assign currents to plane structure
    plane.b = current_span;
    plane.AR = current_AR;
    plane.massPayloads2 = current_cabinetMass;
    plane.nPayloads3 = 3;

    % Wing params change accordingly
    plane.c = plane.b/plane.AR;
    plane.S = plane.b*plane.c;
    plane.croot = 2*plane.c / (1 + plane.taperRatio);
    plane.ctip = plane.croot*plane.taperRatio;
    
    % Tail params change accordingly
    plane = sizeTail(plane);

    % Simulating competition
    [m2results, m3results, plane] = SimulateCompetition(plane, environment);
    
    % Store M2 and M3 results
    plane.m2results = m2results;
    plane.m3results = m3results;
    
    % Store the plane
    planes{i} = plane;
    TOFL(i) = plane.m2results.TOFL;

    % Estimate how long code will take
    if i/lengthCombined >= 0.1 && usedTimer == 0
        estimateTime = msgbox(sprintf('The code will take ~%.2f minutes to run', minutes(seconds(toc) * lengthCombined/i)));
        beep
        usedTimer = 1;
    end
end

%% Calculate scores
[m2scores,m3scores,totalscores,planes] = calcScores(planes);

% Identify best plane
[~, bestPlane_index] = max(totalscores);
bestPlane = planes(bestPlane_index);

% Print best plane
clc; close all; % get rid of loading info
printPlane(bestPlane{1})
beep % just letting the user know the code is done

%% Plotting

% 4D plot about score as a function of span, radar mass, and AR with scores of 0 removed
figure(3)
hold on
grid on
nonZeroTotalScores = totalscores;
nonZeroTotalScores(nonZeroTotalScores == 0) = NaN;
nonZeroTotalScores(nonZeroTotalScores == 3) = NaN;

% scatter3(span(:), AR(:), cabinetMass(:), [], nonZeroTotalScores, 'filled')
% scatter3(span(bestPlane_index), AR(bestPlane_index), cabinetMass(bestPlane_index), [], nonZeroTotalScores(bestPlane_index), 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 2)
% xlabel('Span [m]')
% ylabel('Aspect Ratio [-]')
% zlabel('Cabinet Mass [kg]')
% cb = colorbar;
% if sum(diff(totalscores)) ~= 0
%     clim([0.8*max(totalscores), max(totalscores)])
% end
% cb.Label.String = 'Score [-]';
% view(-37.5, 30)

% Plotting Velocity M2
plotPerformance(bestPlane{1}, 2, 'V', false, '-')
plotPerformance(bestPlane{1}, 2, 'Vstall', true, '--')
plot([0, bestPlane{1}.m2results.performance.t(end)], [bestPlane{1}.Vmax2, bestPlane{1}.Vmax2], 'LineWidth', 2, 'LineStyle', ':')
plotPerformance(bestPlane{1}, 2, 'throttle', true, ':', true, 1)
ylim([0, 1.1])
legend('V', 'V_{stall}', 'V_{target}', 'throttle', 'location', 'southeast')

% Plotting Velocity M3
plotPerformance(bestPlane{1}, 3, 'V', false, '-')
plotPerformance(bestPlane{1}, 3, 'Vstall', true, '--')
plot([0, bestPlane{1}.m3results.performance.t(end)], [bestPlane{1}.Vmax3, bestPlane{1}.Vmax3], 'LineWidth', 2, 'LineStyle', ':')
plotPerformance(bestPlane{1}, 3, 'throttle', true, ':', true, 1)
ylim([0, 1.1])
legend('V', 'V_{stall}', 'V_{target}', 'throttle', 'locaiton', 'southeast')

