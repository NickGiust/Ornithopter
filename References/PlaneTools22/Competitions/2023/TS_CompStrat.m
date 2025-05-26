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

%% Define Plane
plane = definePlane();
plane.AR = 5.5;
plane.b = 1.25; 
%% Define Environment
environment = defineEnvironment("Tucson");

windSpeeds = 0:5:20;
bestPlanes = [];

for w = 1:length(windSpeeds)

%% 5D Trade Study

%Single plane
m2payloadfractionList = 0.30:0.05:0.5;
antennaLengthList = 0.7:0.1:0.9;
GMMultiplierList = 80;

[m2payloadfractionArray,antennaLengthArray,GMMultiplierArray] = ...
    ndgrid(m2payloadfractionList,antennaLengthList, GMMultiplierList);

m2scoreArray = m2payloadfractionArray;
m3scoreArray = m2payloadfractionArray;
planes = {};
boxes = {};
m2scores_USC = m2payloadfractionArray;
m3scores_USC = m2payloadfractionArray;
GMscores_USC = m2payloadfractionArray;
simulationStats = [];
simulationStats.M2TOFLFails = 0;
simulationStats.M3TOFLFails = 0;
simulationStats.invalidTubes = 0;
simulationStats.boxFails = 0;

counter = 0;
totalT = 0;
totalCombos = length(m2payloadfractionList)*length(antennaLengthList)*length(GMMultiplierList);

tic
for p = 1:length(m2payloadfractionList)
    plane.payload2Fraction = m2payloadfractionList(p);

        for l = 1:length(antennaLengthList)
            plane.lengthPayload3 = antennaLengthList(l);
         
                for g = 1:length(GMMultiplierList)

                    plane.m2results = [];
                    plane.m3results = [];
                    plane.GMresults = [];
                    plane.GMloadingMultiplier = GMMultiplierList(g);

                    centerFairingSpan = 0.05;
                    % width of center section (foam?) fairing/wing
                    % section between fuselage and outer wing, m
                    semispan = plane.b/2;
                    outerSectionSpan = semispan - centerFairingSpan - 0.5*plane.wFuse;
                    outerSectionSpanMinimumFraction = 0.83
                    if outerSectionSpan < outerSectionSpanMinimumFraction * semispan
                        outerSectionSpan = outerSectionSpanMinimumFraction * semispan
                    end
                    plane.outerSectionSpan = outerSectionSpan;

                    [plane, box] = checkBoxSizing(plane);

                    plane.validTubes = true;
                    [m2results, m3results, GMresults, plane] = TS_SimulateCompetition(plane, environment, planeToolsDirectory);
                    
                    plane.m2results = m2results;
                    plane.m3results = m3results;
                    plane.GMresults = GMresults;
                    GMMultiplierArray(p,l,g) = plane.GMloadingMultiplier;
                    
                    planes{p,l,g} = plane;
                    boxes{p,l,g} = box;
                    m2scores_USC(p,l,g) = plane.m2results.score;
                    m3scores_USC(p,l,g) = plane.m3results.score;
                    GMscores_USC(p,l,g) = plane.GMresults.score;
                    
                    % Simulation Stats Tracking

                    if ~plane.planeFits_boolean || ~plane.antennaFits_boolean
                        simulationStats.boxFails = simulationStats.boxFails + 1;
                    end

                    if plane.validTubes == false
                        simulationStats.invalidTubes = simulationStats.invalidTubes + 1;
                    end

                    if m2results.TOFLbool == "FAILED"
                        simulationStats.M2TOFLFails = simulationStats.M2TOFLFails + 1;
                    elseif m2results.TOFLbool == "PASSED" && m3results.TOFLbool == "FAILED"
                        simulationStats.M3TOFLFails = simulationStats.M3TOFLFails + 1;
                    end


                    % Simulation Progress Tracking
                    totalT = toc;
                    counter = counter + 1;
                    averageTimePerPlane = totalT/counter;
                    predictedTimeRemaining = averageTimePerPlane*(totalCombos-counter); 
                    
                    fprintf("Total # of Planes in Trade Study: %i\n", totalCombos)
                    fprintf("Planes Completed: %i/%i\n", counter, totalCombos)
                    fprintf("Avg. time per plane: %.2f seconds\n", averageTimePerPlane)
                    fprintf("Estimated time remaining: %.2f minutes\n", predictedTimeRemaining/60)
                    fprintf("Simulation %.2f%% Completed...\n", (counter/totalCombos)*100)
                    
                end
        end
end


%% Best Plane Assumptions, Print Results

totalT

%diary('14x12_6S_14x12_10S.txt');

fprintf("\n--Simulation Stats--\n")
fprintf("Total # of Planes Simulated: %i\n", totalCombos)
fprintf("# of Planes Failing Box Fit: %i\n", simulationStats.boxFails)
fprintf("# of Planes Failing Tube Geometry: %i\n", simulationStats.invalidTubes)
fprintf("# of Planes Failing M2 TOFL: %i\n", simulationStats.M2TOFLFails)
fprintf("# of Planes Failing M3 TOFL: %i\n", simulationStats.M3TOFLFails)
totalFailures = simulationStats.boxFails + simulationStats.invalidTubes + ...
    simulationStats.M2TOFLFails + simulationStats.M3TOFLFails;
fprintf("# of Planes Failing in TOTAL: %i\n", totalFailures)
if totalFailures > totalCombos
    fprintf('\n\n------------------\n')
    fprintf('MAJOR ERROR: MORE PLANES FAILING THAN IN SIMULATION\n')
    fprintf('------------------\n\n\n')
end

% Apply score bonuses (each mission score > 1 possible)
bonuses = true;

% BEST PLANE M2
best_mPayload2 = 8;
best_nLaps2 = 15;
best_m2score = best_mPayload2 * best_nLaps2;

m2scores = m2scores_USC/best_m2score;
% m2scores(m2scores > 1) = 1;

if max(m2scores, [], 'all') > 1
    if bonuses
        m2scores(m2scores > 1) = 1 + (1 - (1./m2scores(m2scores > 1))); % M2 Bonus
    else
        m2scores(m2scores > 1) = 1;
    end
end

% Completion constant
% m2scoresNotZero = m2scores ~= 0;
% m2scores(m2scoresNotZero) = m2scores(m2scoresNotZero) + 1;


% BEST PLANE M3
best_lengthPayload3 = 1.2;
best_time3 = 45;
best_m3score = best_lengthPayload3 / best_time3;

m3scores = m3scores_USC/best_m3score;
%m3scores(m3scores > 1) = 1;

if max(m3scores, [], 'all') > 1
    if bonuses
        m3scores(m3scores > 1) = 1 + (1 - (1./m3scores(m3scores > 1))); % M3 Bonus
    else
        m3scores(m3scores > 1) = 1;
    end
end
% Completion constant
% m3scoresNotZero = m3scores ~= 0;
% m3scores(m2scoresNotZero) = m3scores(m3scoresNotZero) + 2;


% BEST PLANE GM

% best_GMloadingList_lowerLimit = 20;
% best_GMloadingList_upperLimit = 120;
% best_GMloadingList_increment = 5;
% best_GMloadingList = best_GMloadingList_lowerLimit:best_GMloadingList_increment:best_GMloadingList_upperLimit; % SET TO ONE ELEMENT FOR EASY COMPARISON
% bestScoreIndices = best_GMloadingList;

best_GMloadingList = 50;
for i = 1:length(best_GMloadingList)
    
    best_GMloading = best_GMloadingList(i);
    
    GMscores = GMscores_USC/best_GMloading;

    if max(GMscores, [], 'all') > 1
        if bonuses
            GMscores(GMscores > 1) = 1 + (1 - (1./GMscores(GMscores > 1))); % Bonus for exceeding best GM multiplier
        else
            GMscores(GMscores > 1) = 1;
        end
    end
    scoreArray = m2scores + m3scores + GMscores;
    
    % BEST OVERALL SCORING PLANE
    [~, bestScoreIndices(i)] = max(scoreArray(:));

end


bestPlaneIndices = unique(bestScoreIndices, 'stable');

for i = 1:length(bestPlaneIndices)

    % indices of unique plane being best for range of GMbest

    GMIndices = find(bestScoreIndices == bestPlaneIndices(i));
    GMmaxIndex = GMIndices(end);
    GMminIndex = GMIndices(1);

    GMmax = best_GMloadingList(GMmaxIndex);
    GMmin = best_GMloadingList(GMminIndex);

    bestScoreIndex = bestScoreIndices(GMmaxIndex);

    best_plane = planes{bestScoreIndex};
    best_box = boxes{bestScoreIndex};

    m2score = m2scores(bestScoreIndex);
    m3score = m3scores(bestScoreIndex);
    GMscore = GMscores_USC(bestScoreIndex)/GMmax;

    if GMscore > 1
        if bonuses
            GMscore = 1 + (1 - (1/GMscore)); % Bonus for exceeding best GM multiplier
        else
            GMscore = 1;
        end
    end
    
    totalscore = m2score + m3score + GMscore;

    % print results
    
    fprintf("\nFor %.d <= GMbest <= %.d:\n", GMmin, GMmax)
    fprintf("--BEST COMPETITION PLANE--\n")
    
    fprintf("\n--SCORES--\n")
    fprintf("Total Score: %.2f\n", totalscore);
    fprintf("   GM Score: %.2f\n", GMscore);
    fprintf("   M2 Score: %.2f\n", m2score); 
    fprintf("   M3 Score: %.2f\n", m3score); 
    
    fprintf("\n--TAKEOFF (MAX 60 FT)--\n")
    fprintf("M2 TOFL: %s, %.2f ft\n", best_plane.m2results.TOFLbool, best_plane.m2results.TOFL/0.3048);
    TOsegmentIndex = find(best_plane.m2results.performance.V > best_plane.m2results.VTO, 1);
    fprintf("M2 Avg. Current, TO: %.2f A for %.2f s\n", mean(best_plane.m2results.performance.I(1:TOsegmentIndex-1)), best_plane.m2results.coursePoints.startTimes(2));
    fprintf("M3 TOFL: %s, %.2f ft\n", best_plane.m3results.TOFLbool, best_plane.m3results.TOFL/0.3048);
    TOsegmentIndex = find(best_plane.m3results.performance.V > best_plane.m3results.VTO, 1);
    fprintf("M3 Avg. Current, TO: %.2f A for %.2f s\n\n", mean(best_plane.m3results.performance.I(1:TOsegmentIndex-1)), best_plane.m3results.coursePoints.startTimes(2));

    if best_plane.m2results.TOFLbool == "PASSED"
    fprintf("\n--LAPS/PERFORMANCE--\n")
    fprintf("M2: %.d laps in %i minutes %.1f seconds\n    %.1f%% battery energy used\n", best_plane.m2results.nLaps, floor(best_plane.m2results.missionTime/60),rem(best_plane.m2results.missionTime,60),100-100*best_plane.m2results.performance.s(end));
    end
    if best_plane.m3results.TOFLbool == "PASSED"
    fprintf("M3: 3 laps in %.1f seconds\n    %.1f%% battery energy used\n", best_plane.m3results.missionTime, 100-100*best_plane.m3results.performance.s(end));
    end
    fprintf("M2 Max. Cruise Speed: %.0f mph\n", max(best_plane.m2results.performance.V)*2.237)
    fprintf("M3 Max. Cruise Speed: %.0f mph\n", max(best_plane.m3results.performance.V)*2.237)

    fprintf("\n--PAYLOAD--\n")
    fprintf("M2 Payload Fraction: %.2f\n", best_plane.payload2Fraction);
    fprintf("M2 Payload Mass: %.2f kg\n", best_plane.payload2Fraction*best_plane.m2);
    fprintf("M3 Antenna Length: %.2f m\n", best_plane.lengthPayload3);
    fprintf("GM Loading Multiplier (Design Point): %.2f\n", best_plane.GMloadingMultiplier)
    fprintf("GM Loading Multiplier (True): %.2f\n", best_plane.trueGMloadingMultiplier)
    fprintf("GM Load: %.f kg\n", best_plane.trueGMloadingMultiplier*best_plane.m2)
    
    fprintf("\n--MASS--\n")
    fprintf("Empty Mass: %.2f kg\n", best_plane.m_empty);
    fprintf("M2 Mass: %.2f kg\n", best_plane.m2);
    fprintf("M3 Mass: %.2f kg\n", best_plane.m3);
    fprintf("Wing Mass: %.2f kg\n", best_plane.mWing);
    
    fprintf("\n--WING--\n")
    fprintf("Span: %.2f m\n", best_plane.b)
    fprintf("Chord: %.2f m\n", best_plane.c)
    fprintf("AR: %.1f\n", best_plane.AR)

    fprintf("\n--TUBES--\n")
    fprintf("Center Tube Span: %.2f m\n", best_plane.centerTubeSpan)
    fprintf("Center Tube Span Fraction: %.f%%\n", best_plane.centerTubeSpanFraction*100)
    fprintf("Center Tube Outer Diameter: %.2f cm (%.2f in)\n", best_plane.centerTubeOuterDiameter*100, best_plane.centerTubeOuterDiameter/0.0254)
    fprintf("Center Tube Inner Diameter: %.2f cm (%.2f in)\n", best_plane.centerTubeInnerDiameter*100, best_plane.centerTubeInnerDiameter/0.0254)
    fprintf("Center Tube Thickness: %.2f cm (%.2f in)\n", best_plane.centerTubeThickness*100, best_plane.centerTubeThickness/0.0254)
    fprintf("Center Tube Total Carbon Layers: %.d \n\n", best_plane.centerTubeLayers)

    fprintf("Outer Tube Span: %.2f m\n", best_plane.outerSectionSpan)
    fprintf("Outer Tube Min. Inner Diameter: %.2f cm (%.2f in)\n", min(best_plane.outerTubeInnerDiameter)*100, max(best_plane.outerTubeInnerDiameter)/0.0254)
    fprintf("Outer Tube Max Thickness: %.2f cm (%.3f in)\n", max(best_plane.outerTubeThickness)*100, max(best_plane.outerTubeThickness)/0.0254)
    fprintf("Outer Tube Max Carbon Layers: %.d \n", max(best_plane.outerTubeLayers))
    
    fprintf("\n--TAIL--\n")
    fprintf("Span: %.2f m\n", best_plane.bTail)
    fprintf("Chord: %.2f m\n", best_plane.cTail)
    fprintf("Height: %.2f m\n", best_plane.hTail)
    
    fprintf("\n--BOX--\n")
    fprintf("Box length: %.2f in\n", best_box.length/0.0254)
    fprintf("Box width: %.2f in\n", best_box.width/0.0254)
    fprintf("Box height: %.2f in\n", best_box.height/0.0254)
    fprintf("Sum of Dimensions = %.2f in\n\n\n", (best_box.length+best_box.width+best_box.height)/0.0254);

end

bestPlanes = [bestPlanes; best_plane];
end





%{
%% Mass Breakdown
figure
labels = {'Fuselage', 'Wing', 'Tail', 'Battery', 'Motor', 'Landing Gear', 'Servos', 'Wiring', 'ESC'};
pie(best_plane.m_emptyList/best_plane.m_empty, labels);
title("Empty Weight Breakdown")
% ax = gca;
% ax.FontName = "Times New Roman";
% ax.FontSize = 18;
% lgd = legend(labels);
% lgd.Location = 'eastoutside';

%% 2-D PLOTTING

spansAcrossGMbest = [];
ARsAcrossGMbest = [];
USCGMAcrossGMbest = [];
antennaLengthAcrossGMbest = [];
M2payloadFracAcrossGMbest = [];
M2scoresAcrossGMbest = [];
M3scoresAcrossGMbest = [];
GMscoresUSCAcrossGMbest = [];

for i = 1:length(bestScoreIndices)
    spansAcrossGMbest = [spansAcrossGMbest planes{bestScoreIndices(i)}.b];
    ARsAcrossGMbest = [ARsAcrossGMbest planes{bestScoreIndices(i)}.AR];
    M2payloadFracAcrossGMbest = [M2payloadFracAcrossGMbest planes{bestScoreIndices(i)}.payload2Fraction];
    antennaLengthAcrossGMbest = [antennaLengthAcrossGMbest planes{bestScoreIndices(i)}.lengthPayload3];
    USCGMAcrossGMbest = [USCGMAcrossGMbest planes{bestScoreIndices(i)}.trueGMloadingMultiplier];
    M2scoresAcrossGMbest = [M2scoresAcrossGMbest m2scores(bestScoreIndices(i))];
    M3scoresAcrossGMbest = [M3scoresAcrossGMbest m3scores(bestScoreIndices(i))];
    GMscoresUSCAcrossGMbest = [GMscoresUSCAcrossGMbest GMscores_USC(bestScoreIndices(i))];
end

GMscoresAcrossGMbest = GMscoresUSCAcrossGMbest ./ best_GMloadingList;

if max(GMscoresAcrossGMbest, [], 'all') > 1
    if bonuses
        GMscoresAcrossGMbest(GMscoresAcrossGMbest > 1) = ...
            1 + (1 - (1./GMscoresAcrossGMbest(GMscoresAcrossGMbest > 1))); % Bonus for exceeding best GM multiplier
    else
        GMscoresAcrossGMbest(GMscoresAcrossGMbest > 1) = 1;
    end
end

LW = 2; % linewidth for plotting
figure()
subplot(3,1,1)
plot(best_GMloadingList, spansAcrossGMbest, 'LineWidth', LW)
hold on
plot(best_GMloadingList, antennaLengthAcrossGMbest, '-.', 'LineWidth', LW)
plot(best_GMloadingList, M2payloadFracAcrossGMbest, '--', 'LineWidth', LW)
ylim([min(M2payloadFracAcrossGMbest)-0.1 max(spansAcrossGMbest)+0.1])
ax = gca;
ax.YAxis.TickLabelFormat = "%.1f";
legend('b [m]', 'L_{antenna} [m]', 'M2 Payload Fraction [-]', 'Location', 'best')
xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
grid on
ax = gca;
ax.FontSize = 14;

subplot(3,1,2)
yyaxis left
plot(best_GMloadingList, ARsAcrossGMbest, 'LineWidth', LW)
ylim([min(ARsAcrossGMbest) - 1, max(ARsAcrossGMbest) + 1])
xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
ylabel('AR [-]')
ax = gca;

yyaxis right
plot(best_GMloadingList, USCGMAcrossGMbest, '*-.', 'LineWidth', LW)
ylabel('USC GM Margin [-]')
ax = gca;
ax.FontSize = 14;



subplot(3,1,3)
plot(best_GMloadingList, M3scoresAcrossGMbest, '-.', 'LineWidth', LW)
hold on
plot(best_GMloadingList, M2scoresAcrossGMbest, '--', 'LineWidth', LW)
plot(best_GMloadingList, GMscoresAcrossGMbest, '-', 'LineWidth', LW)

% ylim([0 1])
ax = gca;
ax.YAxis.TickLabelFormat = "%.1f";
ax.FontSize = 14;
legend('M3 Score', 'M2 Score', 'GM Score', 'location', 'best')
xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
xlabel('Best Team GM Margin [-]')
grid on


%% 3-D PLOTTING

% close all;

hiddenLengthIndices = [];
for p = 1:length(m2payloadfractionList)
    for s = 1:length(spanList)
        for a = 1:length(ARList)
            for g = 1:length(GMMultiplierList)
                [M,I] = max(scoreArray(p,s,:,a,g));
                hiddenLengthIndices = [hiddenLengthIndices I];
            end
        end
    end
end


for g = 1:length(GMMultiplierList)
    figure
    hold on
    for p = 1:length(m2payloadfractionList)
        for i = 1:hiddenLengthIndices
            span = spanArray(p,:,hiddenLengthIndices(i),:,g);
            span = squeeze(span);
            AR = ARArray(p,:,hiddenLengthIndices(i),:,g);
            AR = squeeze(AR);
%             GMload = GMMultiplierArray(p,:,hiddenLengthIndices(i),:,g);
%             GMload = squeeze(GMload);
            m2payloadFraction = m2payloadfractionArray(p,:,hiddenLengthIndices(i),:,g);
            m2payloadFraction = squeeze(m2payloadFraction);
            score = scoreArray(p,:,hiddenLengthIndices(i),:,g);
            score = squeeze(score);
            %scatter3(span(:), AR(:), GMload(:), 80, score(:), 'filled', 'MarkerFaceAlpha', 0.2);
            surf(span, AR, m2payloadFraction, score)
        end
    end
    view(3)
    % shading interp
    alpha 0.5
    titleString = strcat("USC Design GM Load Margin: ", num2str(GMMultiplierList(g)));
    title(titleString,'FontSize',18,'FontWeight','bold')
    xlabel("Span [m]",'FontSize',18,'FontWeight','bold')
    ylabel("AR [-]",'FontSize',18,'FontWeight','bold')
    zlabel("M2 Payload Fraction [-]",'FontSize',18,'FontWeight','bold')
    cb = colorbar; % create and label the colorbar
    cb.Label.String = 'Competition Score';
    cb.Label.FontSize = 18;
    cb.Label.FontWeight = 'bold';
    cb.Limits = ([0 2.5]);
    grid on
    lighting phong
%         set('edgecolor',[0 0 0.4],'meshstyle','both','linewidth',.15);

end
%}

