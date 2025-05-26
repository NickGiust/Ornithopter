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
environment = defineEnvironment("Whittier");

%% Trade Study Setup

% Define parameter ranges and values
span_range = 1.18:0.1:1.68;        % [m] max is 1.68 m due to material constraints
AR_range = 1:0.25:7;              % [-]
radarMass_range = 1:0.25:2.75; % [kg]

% span_range = 1.68;    % [m] max is 1.68 m due to material constraints
% AR_range = 4.5;        % [-]
% radarMass_range = 0.5; % [kg]

% Generate grids of parameter values
[span, AR, radarMass] = meshgrid(span_range, AR_range, radarMass_range);

% Generate grids of indices
[span_index, AR_index, radarMass_index] = ndgrid(1:numel(span_range), 1:numel(AR_range), 1:numel(radarMass_range));

% Combine parameter grids and indices
combined = cat(2, reshape(span, [], 1), reshape(AR, [], 1), reshape(radarMass, [], 1), ...
                   reshape(span_index, [], 1), reshape(AR_index, [], 1), reshape(radarMass_index, [], 1));
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
    current_radarMass = combined(i, 3);
%     current_span_index = combined(i, 4);        % not used
%     current_AR_index2 = combined(i, 5);         % not used
%     current_radarMass_index3 = combined(i, 6);  % not used

    % Print the current plane specs so user knows it's running
    fprintf('Plane:            %.0f of %.0f\nRadar mass [kg]:  %.2f\nSpan [m]:         %.1f\nAR:               %.1f\n\n',...
        i, lengthCombined, current_radarMass, current_span, current_AR);
    
    % Assign currents to plane structure
    plane.b = current_span;
    plane.AR = current_AR;
    plane.radarMass = current_radarMass;

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

% % 4D plot about score as a function of span, radar mass, and AR (includes failed TOFL flights)
% figure(2)
% scatter3(span(:), AR(:), radarMass(:), [], totalscores(:), 'filled')
% xlabel('Span [m]')
% ylabel('Aspect Ratio [-]')
% zlabel('Radar Mass [kg]')
% cb = colorbar;
% if sum(diff(totalscores)) ~= 0
%     clim([0.8*max(totalscores), max(totalscores)])
% end
% % cb.Limits = [0.8*max(totalscores), max(totalscores)];
% colormap(cb, parula);
% cb.Label.String = 'Overall Competition Score';

% 4D plot about score as a function of span, radar mass, and AR with scores of 0 removed
figure(3)
hold on
grid on
nonZeroTotalscores = totalscores;
nonZeroTotalscores(nonZeroTotalscores == 0) = NaN;

scatter3(span(:), AR(:), radarMass(:), [], nonZeroTotalscores(:), 'filled')
scatter3(span(bestPlane_index), AR(bestPlane_index), radarMass(bestPlane_index), [], nonZeroTotalscores(bestPlane_index), 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 2)
xlabel('Span [m]')
ylabel('Aspect Ratio [-]')
zlabel('Radar Mass [kg]')
cb = colorbar;
if sum(diff(totalscores)) ~= 0
    clim([0.8*max(totalscores), max(totalscores)])
end
cb.Label.String = 'Overall Competition Score';
view(-37.5, 30)
