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
environment = defineEnvironment("Wichita");

%% Trade Study Setup

% Define parameter ranges and values
% overall
span_range = 1.5;  % [m] max is 1.52 m (5ft)
AR_range = 1.5:0.25:5;             % [-]
% M2
cabinetMass_range = 0.5:0.1:2.5; % [kg]
% M3
nPassengers_range = 15:5:40; % [-] only whole numbers

% Define parameter ranges and values
% overall
span_range = 1.5;  % [m] max is 1.52 m (5ft)
AR_range = 3.8;             % [-]
% M2
cabinetMass_range = 1.4; % [kg]
% M3
nPassengers_range = 35; % [-] only whole numbers

% % for testing prop packages
% AR_range = 3.5;
% cabinetMass_range = 1.5;
% nPassengers_range = 40;
% motorPropVal = {'TMotor_AT4125_540Kv',           14, 8, 1;...
%                 'Scorpion_HKIII_5025_520KV_F3S', 14, 8, 1;...
%                 'Hacker_A50_12L_V3',             15, 8, 1;...
%                 'Scorpion_A_4220_540',           11, 7, 2;...
%                 'Scorpion_A_4220_540',           11, 6, 2;...
%                 'Scorpion_SII_4020_540KV',       11, 7, 2;...
%                 'Scorpion_SII_4020_540KV',       11, 6, 2};
% motorProp_range = 1:size(motorPropVal, 1);
% propellerMAT = matfile('myFile2.mat');
% allProps = propellerMAT.allProps;

% Generate grids of parameter values
[AR, cabinetMass, nPassengers] = ndgrid(AR_range, cabinetMass_range, nPassengers_range);

% Generate grids of indices
[AR_index, cabinetMass_index, nPassengers_index] = ndgrid(1:numel(AR_range), 1:numel(cabinetMass_range), 1:numel(nPassengers_range));

% Combine parameter grids and indices
combined = cat(2, reshape(AR, [], 1), reshape(cabinetMass, [], 1), reshape(nPassengers, [], 1),...
                   reshape(AR_index, [], 1), reshape(cabinetMass_index, [], 1), reshape(nPassengers_index, [], 1));
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
    current_AR = combined(i, 1);
    current_cabinetMass = combined(i, 2);
    current_nPassengers = combined(i, 3);
    current_span = span_range;

%     Print the current plane specs so user knows it's running
    fprintf('Plane:              %.0f of %.0f\nSpan [m]:           %.2f\nCabinet mass [kg]:  %.2f\nNum Passengers [-]: %.2f\nAR:                 %.1f\n\n',...
       i, lengthCombined, current_span, current_cabinetMass, current_nPassengers, current_AR);
    
%     fprintf('Motor:   %s\nDiam:    %.0f\nPitch:   %.0f\nnMotors: %.0f\n\n', motorPropVal{i,1}, motorPropVal{i,2}, motorPropVal{i,3}, motorPropVal{i,4})

    % Assign currents to plane structure
    plane.AR = current_AR;
    plane.massPayloads2 = current_cabinetMass;
    plane.nPayloads3 = current_nPassengers;
    plane.b = current_span;
    
%     % change motor and prop
%     plane.motorType = motorPropVal{i, 1};
%     warning('off')
%     load('motor.mat', plane.motorType);
%     motor = eval(plane.motorType);
%     plane.Kv = motor.Kv;
%     plane.motorMaxPower = motor.maxPower;
%     Rmotor = motor.R;
%     plane.motorLength = motor.length;
%     plane.motorDiam = motor.diam;
% 
%     Resc = .0005;
%     Rwire = .01;
%     plane.Rt2 = Resc + Rwire + Rmotor + (plane.bat2R/plane.nParallel2);
%     plane.Rt3 = Resc + Rwire + Rmotor + (plane.bat3R/plane.nParallel3);
% 
%     plane.nMotors = motorPropVal{i, 4};
% 
%     diameter = motorPropVal{i, 2};
%     pitch = motorPropVal{i, 3};
%     diamIndices = find(allProps(:,1)==diameter);
%     pitchIndices = find(allProps(:,2)==pitch);
%     [index,~]=intersect(diamIndices,pitchIndices);
%     
%     plane.D2 = diameter*0.0254;
%     plane.P2 = pitch;
%     plane.A2 = allProps(index,3);
%     plane.B2 = allProps(index,4);
%     plane.C2 = allProps(index,5);
%     plane.E2 = allProps(index,6);
%     plane.F2 = allProps(index,7);
%     plane.G2 = allProps(index,8);
% 
%     plane.D3 = diameter*0.0254;
%     plane.P3 = pitch;
%     plane.A3 = allProps(index,3);
%     plane.B3 = allProps(index,4);
%     plane.C3 = allProps(index,5);
%     plane.E3 = allProps(index,6);
%     plane.F3 = allProps(index,7);
%     plane.G3 = allProps(index,8);

    % Fuselage width changes according to payload
    % width of 1 patient (1" + 11/16") + 1 emt (1.5") (plus some fudge room) [in -> m]
    wFuse_m2lim = ((1 + 11/16) + 1.5) * 0.0254 * 1.3;
    % (width of pax) * (num Pax across) * (factor of fudge) [in -> m]
    nPaxWide = 6; % number of passengers across the width
    wFuse_m3lim = 1.5 * nPaxWide * 0.0254 * 1.3;
    plane.wFuse = max([wFuse_m2lim, wFuse_m3lim]);

    % Fuselage length changes according to payload
    lFuse_noseCone = 1.2 * plane.hFuse;
    lFuse_tailCone = 3 * plane.hFuse;
    % (width crew + length gourney + length cabinet) * (factor of fudge) [in -> m]
    lFuse_m2lim = (1.5 + 5.5 + 3.0) * 1.3 * 0.0254 + lFuse_noseCone + lFuse_tailCone;
    % (width crew + width pax / X pax per row) * (factor of fudge) [in -> m]
    lFuse_m3lim = (1.5 + ceil(current_nPassengers/nPaxWide) * 1.5) * 1.3 * 0.0254 + lFuse_noseCone + lFuse_tailCone;
    plane.lFuse = max([lFuse_m2lim, lFuse_m3lim]); % pick the larger of the 2 fuselages sizes for the given plane

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
fprintf('--INPUTS--\n')
fprintf('AR: %.1f\n', bestPlane{1}.AR)
fprintf('Cabinet Mass: %.2f\n', bestPlane{1}.massPayloads2)
fprintf('Num Pax: %.0f\n', bestPlane{1}.nPayloads3)
printPlane(bestPlane{1})
beep % just letting the user know the code is done

%% Plotting

% % plotting TOFL as a function of prop package
% figure()
% hold on
% TOFLVal = zeros(length(planes), 2);
% motorName = zeros(1, length(planes));
% massM2 = cell(1, length(planes));
% massM3 = cell(1, length(planes));
% 
% for i = 1:length(planes)
%     TOFLVal(i, :) = [convlength(planes{i}.m2results.TOFL, 'm', 'ft'),...
%                      convlength(planes{i}.m3results.TOFL, 'm', 'ft')];
%     motorName(i) = i;
%     massM2{i} = sprintf('%.2f kg', planes{i}.m2);
%     massM3{i} = sprintf('%.2f kg', planes{i}.m3);
% end
% 
% bGraph = bar(motorName, TOFLVal);
% yline(20, 'LineWidth', 2, 'LineStyle', ':', 'Color', '#d3d3d3')
% 
% % adding masses
% xtipsm2 = bGraph(1).XEndPoints;
% ytipsm2 = bGraph(1).YEndPoints;
% text(xtipsm2,ytipsm2,massM2,'HorizontalAlignment','center','VerticalAlignment','bottom')
% xtipsm3 = bGraph(2).XEndPoints;
% ytipsm3 = bGraph(2).YEndPoints;
% text(xtipsm3,ytipsm3,massM3,'HorizontalAlignment','center','VerticalAlignment','bottom')
% 
% % formatting
% bGraph(1).FaceColor = '#990000';
% bGraph(2).FaceColor = '#FFC72C';
% ylabel('TOFL [ft]')
% title('Trading Prop Packages for TOFL')
% legend('M2', 'M3', 'Location', 'northwest')

% 4D plot about score as a function of span, radar mass, and AR with scores of 0 removed
figure()
hold on
grid on
nonZeroTotalScores = totalscores;
nonZeroTotalScores(nonZeroTotalScores == 0) = NaN; % remove any scores that are 0 (failed M2 and M3)
nonZeroTotalScores_index = ~isnan(nonZeroTotalScores);
M2OnlyScores = nonZeroTotalScores;
M2OnlyScores(nonZeroTotalScores > 3) = NaN;
M2OnlyScores_idx = ~isnan(M2OnlyScores);
M3OnlyScores = nonZeroTotalScores;
M3OnlyScores(nonZeroTotalScores < 3) = NaN; % remove any scores that failed M3 TOFL (score is under 3)

scatter3(nPassengers(:), AR(:), cabinetMass(:), [], M3OnlyScores(:), 'filled')
scatter3(nPassengers(nonZeroTotalScores_index), AR(nonZeroTotalScores_index), cabinetMass(nonZeroTotalScores_index), 'MarkerEdgeColor', 'black')
% scatter3(span(bestPlane_index), AR(bestPlane_index), cabinetMass(bestPlane_index), [], nonZeroTotalScores(bestPlane_index), 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 2)
xlabel('nPassengers [-]')
ylabel('Aspect Ratio [-]')
zlabel('Cabinet Mass [kg]')
cb = colorbar;
if sum(diff(totalscores)) ~= 0
    clim([0.9*max(totalscores), max(totalscores)])
end
cb.Label.String = 'Score [-]';
set(gca,'Ydir','reverse')
colormap('jet');

% 3d plot
view(-37.5, 30)
% zlim([0.5 1.5])
% xlim([1.2 1.5])
% ylim([2 5])

% % slices
% xView = 1.3;
% titleText = sprintf('Span = %.1f', xView);
% title(titleText)
% view(-90, 0)
% xlim([xView - 0.01, xView + 0.01])

% trade study: nPax, AR, mCabinet
% plotting surfaces to show score

figure()
hold on
zlabel('Number of Pax [-]')
xlabel('Aspect Ratio [-]')
ylabel('Cabinet Mass [kg]')
M3OnlyScoresMatrix = reshape(M3OnlyScores, [length(AR_range), length(cabinetMass_range), length(nPassengers_range)]);
M2OnlyScoresIdxMatrix = reshape(M2OnlyScores_idx, [length(AR_range), length(cabinetMass_range), length(nPassengers_range)]);

for i = 1:length(nPassengers_range)
    constPaxIdx = find(nPassengers == nPassengers_range(i));

    surf(AR(:,:,i), cabinetMass(:,:,i), nPassengers(:,:,i), M3OnlyScoresMatrix(:,:,i),...
        'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'interp')

%     surf(AR(:,:,i), cabinetMass(:,:,i), nPassengers(:,:,i), M2ScoreColor(:,:,i),...
%         'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', '#DEDEDE')
end

cb = colorbar;
if sum(diff(totalscores)) ~= 0
    clim([0.9*max(totalscores), max(totalscores)])
end
cb.Label.String = 'Score [-]';
colormap('jet');
view(-75, 20)
ylim([0.5, 1.5])
xlim([2, 5])
zlim([15, 40])
grid on

% adding on the design point for V0
scatter3(3.8, 1.4, 35, 100, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k')
scatter3(bestPlane{1}.AR, bestPlane{1}.massPayloads2, bestPlane{1}.nPayloads3, 100, 'red', 'o', 'LineWidth', 2)

% Plotting Velocity M2
plotPerformance(bestPlane{1}, 2, 'V', false, '-')
plotPerformance(bestPlane{1}, 2, 'Vstall', true, '--')
plot([0, bestPlane{1}.m2results.performance.t(end)], [1.05*bestPlane{1}.Vmax2, 1.05*bestPlane{1}.Vmax2], 'LineWidth', 2, 'LineStyle', ':')
plot([0, bestPlane{1}.m2results.performance.t(end)], [0.95*bestPlane{1}.Vmax2, 0.95*bestPlane{1}.Vmax2], 'LineWidth', 2, 'LineStyle', ':')
plotPerformance(bestPlane{1}, 2, 'throttle', true, ':', true, 1)
ylim([0, 1.1])
legend('V', 'V_{stall}', 'V_{target,max}', 'V_{target,min}', 'throttle', 'location', 'southeast')

% Plotting Velocity M3
plotPerformance(bestPlane{1}, 3, 'V', false, '-')
plotPerformance(bestPlane{1}, 3, 'Vstall', true, '--')
if isfield(bestPlane{1}.m3results, "performance")
    plot([0, bestPlane{1}.m3results.performance.t(end)], [1.05*bestPlane{1}.Vmax3, 1.05*bestPlane{1}.Vmax3], 'LineWidth', 2, 'LineStyle', ':')
    plot([0, bestPlane{1}.m3results.performance.t(end)], [0.95*bestPlane{1}.Vmax3, 0.95*bestPlane{1}.Vmax3], 'LineWidth', 2, 'LineStyle', ':')
end
plotPerformance(bestPlane{1}, 3, 'throttle', true, ':', true, 1)
ylim([0, 1.1])
legend('V', 'V_{stall}', 'V_{target}', 'throttle', 'locaiton', 'southeast')

% Plotting radius of turns M2
plotPerformance(bestPlane{1}, 2, 'r', false, '-')
ylim([0, 30])

plotPerformance(bestPlane{1}, 2, 'throttle')
ylim([0,1.1])

