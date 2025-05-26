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
slash
planeToolsDirectory = [currentDirectory{1} 'PlaneTools22'];

%% Add directories to path and initialize global propulsion functions
addpath([planeToolsDirectory slash 'PhysicsModels' slash]);
addpath([planeToolsDirectory slash 'PhysicsModels' slash 'MassBuildup']);
addpath([planeToolsDirectory slash 'ComponentLibrary']);
addpath([planeToolsDirectory slash 'Plotting']);
addpath([planeToolsDirectory slash 'Environments'])
% addpath([planeToolsDirectory slash 'ComponentLibrary' slash 'material.mat']);


%% 5D Trade Study

environment = defineEnvironment("Tucson");

% m2payloadfractionList = 0.3:0.05:0.7;
% ARList = 4:0.5:10;
% spanList = 0.8:0.05:2.0;
% antennaLengthList = [0.90 0.95];
% GMMultiplierList = 40:10:160;

% SMALL STUDY
% m2payloadfractionList = 0.3:0.05:0.55;
% ARList = 5:0.5:8;
% spanList = 1.2:0.05:1.7;
% antennaLengthList = 0.85;
% GMMultiplierList = 60:10:90;

%Single plane
m2payloadfractionList = 0.43;
ARList = 7.2;
spanList = 1.55; 
antennaLengthList = 0.85;
GMMultiplierList = 80;


m2payloadfractionList = 0.31:0.02:0.65;
spanList = 1.2:0.05:2.1;
ARList = 6; % doesn't quite work if trying to hold c constant
antennaLengthList = 0.85;
GMMultiplierList = 50:5:80;


[m2payloadfractionArray,spanArray,antennaLengthArray,ARArray,GMMultiplierArray] = ...
    ndgrid(m2payloadfractionList,spanList,antennaLengthList,ARList, GMMultiplierList);

m2scoreArray = zeros(size(spanArray));
m3scoreArray = m2scoreArray;
planes = {};
boxes = {};
m2scores_USC = m2scoreArray;
m3scores_USC = m2scoreArray;
GMscores_USC = m2scoreArray;
M2TOFLFails = m2scoreArray;
M3TOFLFails = m2scoreArray;
invalidTubes = m2scoreArray;
boxFails = m2scoreArray;
completedSim = m2scoreArray;

counter = 0;
totalT = 0;
totalCombos = length(m2payloadfractionList)*length(spanList)*length(antennaLengthList)*length(GMMultiplierList)*length(ARList)

length_p = length(m2payloadfractionList);
length_s = length(spanList);
length_l = length(antennaLengthList);
length_a = length(ARList);
length_g = length(GMMultiplierList);

tic
parfor p = 1:length_p

    for s = 1:length_s

        for l = 1:length_l

            for a = 1:length_a

                for g = 1:length_g

                    plane = definePlane();

                    plane.payload2Fraction = m2payloadfractionList(p);
    
                    plane.b = spanList(s);

                    centerFairingSpan = 0.035;
                    % width of center section (foam?) fairing/wing
                    % section between fuselage and outer wing, m
                    semispan = plane.b/2;
                    outerSectionSpan = semispan - centerFairingSpan - 0.5*plane.wFuse;
                    outerSectionSpanMinimumFraction = 0.82
                    if outerSectionSpan < outerSectionSpanMinimumFraction * semispan
                        outerSectionSpan = outerSectionSpanMinimumFraction * semispan
                    end
                    plane.outerSectionSpan = outerSectionSpan;
                        
                    plane.lengthPayload3 = antennaLengthList(l);
    
                    %{
                    plane.AR = ARList(a);
                    plane.c = plane.b/plane.AR;
                    %}
                    plane.c = 0.23;
                    plane.AR = plane.b/plane.c;
                    plane.S = plane.b*plane.c;
                    % MATH
                    % cavg = 0.5*(croot + ctip)
                    % taperRatio = ctip/croot
                    % ctip = croot*taperRatio
                    % cavg = 0.5*(croot + croot*taperRatio)
                    % cavg = 0.5*croot*(1+taperRatio)
                    % croot = 2*cavg / (1 + taperRatio)
                    plane.croot = 2*plane.c / (1 + plane.taperRatio);
                    plane.ctip = plane.croot*plane.taperRatio;
    
                    plane.lh = 0.82*plane.lengthPayload3 % boom is shorter than antenna to ensure box fitment
                    plane.Sh = plane.Vh*plane.S*plane.c/plane.lh;
                    plane.bh = sqrt(plane.ARh*plane.Sh)
                    plane.ch = plane.Sh./plane.bh
                    plane.lv = plane.lh;
                    plane.Sv = plane.Vv*plane.S*plane.b/plane.lv;
                    plane.cv = plane.ch;
                    plane.bv = plane.Sv/plane.cv;
                    plane.ARv = plane.bv.^2./plane.Sv;
                    % Tail
                    plane.cTail = plane.ch;
                    plane.bTail = plane.bh;
                    plane.hTail = plane.bv;



                    plane.m2results = [];
                    plane.m3results = [];
                    plane.GMresults = [];
                    plane.GMloadingMultiplier = GMMultiplierList(g);

                    [plane, box] = checkBoxSizing(plane);

                    plane.validTubes = true;
                    [m2results, m3results, GMresults, plane] = TS_SimulateCompetition(plane, environment, planeToolsDirectory);
                    
                    plane.m2results = m2results;
                    plane.m3results = m3results;
                    plane.GMresults = GMresults;
%                     GMMultiplierArray(p,s,l,a,g) = plane.GMloadingMultiplier;
                    
                    planes{p,s,l,a,g} = plane;
                    boxes{p,s,l,a,g} = box;
                    m2scores_USC(p,s,l,a,g) = plane.m2results.score;
                    m3scores_USC(p,s,l,a,g) = plane.m3results.score;
                    GMscores_USC(p,s,l,a,g) = plane.GMresults.score;

                    
                    % Simulation Stats Tracking

                    
                    if ~plane.planeFits_boolean || ~plane.antennaFits_boolean
                        boxFails(p,s,l,a,g) = 1;
                    end

                    if plane.validTubes == false
                        invalidTubes(p,s,l,a,g) = 1;
                    end

                    if m2results.TOFLbool == "FAILED"
                        M2TOFLFails(p,s,l,a,g) = 1;
                    elseif m2results.TOFLbool == "PASSED" && m3results.TOFLbool == "FAILED"
                        M3TOFLFails(p,s,l,a,g) = 1;
                    end

%                     completedSim(p,s,l,a,g) = 1;
%                     disp('test')
%                     planesCompleted = sum(completedSim);
% 
                    clc;
%                     fprintf("Total # of Planes in Trade Study: %i\n", totalCombos)
%                     fprintf("Planes Completed: %i/%i\n", planesCompleted, totalCombos)
%                     fprintf("Avg. time per plane: %.2f seconds\n", averageTimePerPlane)
%                     fprintf("Estimated time remaining: %.2f minutes\n", predictedTimeRemaining/60)
%                     fprintf("Simulation %.2f%% Completed...\n", (planesCompleted/totalCombos)*100)
%                     

%{
                    % Simulation Progress Tracking
                    totalT = toc;
                    counter = counter + 1;
                    averageTimePerPlane = totalT/counter;
                    predictedTimeRemaining = averageTimePerPlane*(totalCombos-counter); 
                    
                    clc;
                    fprintf("Total # of Planes in Trade Study: %i\n", totalCombos)
                    fprintf("Planes Completed: %i/%i\n", counter, totalCombos)
                    fprintf("Avg. time per plane: %.2f seconds\n", averageTimePerPlane)
                    fprintf("Estimated time remaining: %.2f minutes\n", predictedTimeRemaining/60)
                    fprintf("Simulation %.2f%% Completed...\n", (counter/totalCombos)*100)
%}
                    
                end
            end
        end
    end
end


totalT = toc
%% Best Plane Assumptions, Print Results

% clc

%diary('14x12_6S_14x12_10S.txt'); % use clock 


totalBoxFails = sum(boxFails, 'all');
totalInvalidTubes = sum(invalidTubes, 'all');
totalM2TOFLFails = sum(M2TOFLFails, 'all');
totalM3TOFLFails = sum(M3TOFLFails, 'all');
totalFailures = totalBoxFails + totalInvalidTubes + ...
    totalM2TOFLFails + totalM3TOFLFails;

fprintf("\n--Simulation Stats--\n")
fprintf("Total # of Planes Simulated: %i\n", totalCombos)
fprintf("# of Planes Failing Box Fit: %i\n", totalBoxFails)
fprintf("# of Planes Failing Tube Geometry: %i\n", totalInvalidTubes)
fprintf("# of Planes Failing M2 TOFL: %i\n", totalM2TOFLFails)
fprintf("# of Planes Failing M3 TOFL: %i\n", totalM3TOFLFails)

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
best_nLaps2 = 16;
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
best_lengthPayload3 = 1.1;
best_time3 = 50;
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

best_GMloadingList_lowerLimit = 50;
best_GMloadingList_upperLimit = 50;
best_GMloadingList_increment = 1;
best_GMloadingList = best_GMloadingList_lowerLimit:best_GMloadingList_increment:best_GMloadingList_upperLimit; % SET TO ONE ELEMENT FOR EASY COMPARISON
bestScoreIndices = best_GMloadingList;

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
    fprintf("M2 TO: %s, %.2f ft TOFL\n", best_plane.m2results.TOFLbool, best_plane.m2results.TOFL/0.3048);
    if isfield(best_plane.m2results, 'performance')
        TOsegmentIndex = find(best_plane.m2results.performance.V > best_plane.m2results.VTO, 1);
        fprintf("       Avg. Current, TO: %.2f A for %.2f s\n", mean(best_plane.m2results.performance.I(1:TOsegmentIndex-1)), best_plane.m2results.coursePoints.startTimes(2));
        fprintf("       VStall: %.2f m/s, VTO: %.2f m/s\n", best_plane.m2results.Vstall, best_plane.m2results.VTO)
        fprintf("M3 TO: %s, %.2f ft TOFL\n", best_plane.m3results.TOFLbool, best_plane.m3results.TOFL/0.3048);
        if isfield(best_plane.m3results, 'performance')
            TOsegmentIndex = find(best_plane.m3results.performance.V > best_plane.m3results.VTO, 1);
            fprintf("       Avg. Current, TO: %.2f A for %.2f s\n", mean(best_plane.m3results.performance.I(1:TOsegmentIndex-1)), best_plane.m3results.coursePoints.startTimes(2));
            fprintf("       Stall: %.2f m/s, VTO: %.2f m/s\n\n", best_plane.m3results.Vstall, best_plane.m3results.VTO)
        end
    end

    if best_plane.m2results.TOFLbool == "PASSED"
        fprintf("\n--LAPS/PERFORMANCE--\n")
        fprintf("M2: %.d laps in %i minutes %.1f seconds\n    %.1f%% battery energy used\n", best_plane.m2results.nLaps, floor(best_plane.m2results.missionTime/60),rem(best_plane.m2results.missionTime,60),100-100*best_plane.m2results.performance.s(end));
        fprintf("    Max. Cruise Speed: %.0f mph\n", max(best_plane.m2results.performance.V)*2.237)
    end
    if best_plane.m3results.TOFLbool == "PASSED"
        fprintf("M3: 3 laps in %.1f seconds\n    %.1f%% battery energy used\n", best_plane.m3results.missionTime, 100-100*best_plane.m3results.performance.s(end));
        fprintf("    Max. Cruise Speed: %.0f mph\n", max(best_plane.m3results.performance.V)*2.237)
    end
    

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
    fprintf("Chord: %.3f m\n", best_plane.c)
    fprintf("AR: %.1f\n", best_plane.AR)
    fprintf("S: %.3f m^2\n", best_plane.S)

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
    fprintf("Span: %.3f m\n", best_plane.bTail)
    fprintf("Chord: %.3f m\n", best_plane.cTail)
    fprintf("Height: %.3f m\n", best_plane.hTail)
    
    fprintf("\n--BOX--\n")
    fprintf("Box is %s design\n", best_box.type)
    fprintf("Box length: %.2f in\n", best_box.length/0.0254)
    fprintf("Box width: %.2f in\n", best_box.width/0.0254)
    fprintf("Box height: %.2f in\n", best_box.height/0.0254)
    fprintf("Sum of Dimensions = %.2f in\n\n\n", (best_box.length+best_box.width+best_box.height)/0.0254);

end



%% Drag breakdown pie charts

figure()
t = tiledlayout(1,2,'TileSpacing','compact');

% Create pie charts
ax1 = nexttile;
pie(ax1,best_plane.dragResults2.DAshares*100)
title('Mission 2')
ax = gca;
ax.FontSize = 14;

ax2 = nexttile;
pie(ax2,best_plane.dragResults3.DAshares*100)
title('Mission 3')
ax = gca;
ax.FontSize = 14;

% Create legend
lgd = legend(["Wing", "", "Fuselage", "Tail", "Antenna", "Landing Gear"]);
lgd.Layout.Tile = 'east';


%{
figure()
pie(best_plane.dragResults2.DAshares*100, ["Wing", "", "Fuselage", "Tail", "Payloads", "Landing Gear"]);
ax = gca;
ax.FontSize = 14;

figure()
pie(best_plane.dragResults3.DAshares*100, ["Wing", "", "Fuselage", "Tail", "Antenna", "Landing Gear"]);
ax = gca;
ax.FontSize = 14;
%}


%% 2-D PLOTTING



% close all
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

% Chose unit system to plot in here!
% units = "metric";
units = "metric";


if units == "metric"

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
    plot(best_GMloadingList, USCGMAcrossGMbest, '-.', 'LineWidth', LW)
    ylabel('USC GM Load Margin [-]')
    ax = gca;
    ax.FontSize = 14;
    
    
    subplot(3,1,3)
    plot(best_GMloadingList, GMscoresAcrossGMbest, '-', 'LineWidth', LW)
    hold on
    plot(best_GMloadingList, M3scoresAcrossGMbest, '-.', 'LineWidth', LW)
    plot(best_GMloadingList, M2scoresAcrossGMbest, '--', 'LineWidth', LW)
    
    ax = gca;
    ax.YAxis.TickLabelFormat = "%.1f";
    ax.FontSize = 14;
    ylim([0 2])
    legend('GM Score', 'M3 Score', 'M2 Score', 'location', 'best')
    xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
    xlabel('Best Team GM Load Margin [-]')
    grid on
    
    % Mass Breakdown
    figure
    labels = {'Fuselage', 'Wing', 'Tail', 'Battery', 'Motor', 'Landing Gear', 'Servos', 'Wiring', 'ESC'};
    pie(best_plane.m_emptyList/best_plane.m_empty, labels);
    title("Empty Weight Breakdown")
    % ax = gca;
    % ax.FontName = "Times New Roman";
    % ax.FontSize = 18;
    % lgd = legend(labels);
    % lgd.Location = 'eastoutside';

elseif units == "imperial"

    figure()
    subplot(3,1,1)
    plot(best_GMloadingList, ARsAcrossGMbest, 'LineWidth', LW)    
    hold on
    plot(best_GMloadingList, spansAcrossGMbest/0.3048, '-.', 'LineWidth', LW)
    ylim([min([spansAcrossGMbest/0.3048 ARsAcrossGMbest])-0.5 max([spansAcrossGMbest/0.3048 ARsAcrossGMbest])+0.5])
    ax = gca;
    ax.YAxis.TickLabelFormat = "%.f";
    legend('AR [-]', 'b [ft]', 'Location', 'best')
    xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
    grid on
    ax = gca;
    ax.FontSize = 14;
    
    subplot(3,1,2)
    yyaxis left
    plot(best_GMloadingList, antennaLengthAcrossGMbest/0.3048, 'LineWidth', LW)
    hold on
    plot(best_GMloadingList, M2payloadFracAcrossGMbest, '--', 'LineWidth', LW)
    ylim([0 max(antennaLengthAcrossGMbest/0.3048) + 0.5])
    xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
    grid on
    
    yyaxis right
    plot(best_GMloadingList, USCGMAcrossGMbest, '-.', 'LineWidth', LW)
    ylabel('USC GM Load Margin [-]')
    ax = gca;
    ax.FontSize = 14;
    legend('L_{antenna} [ft]', 'M2 Payload Fraction [-]', 'USC GM Load Margin', 'Location', 'best')

    
    
    subplot(3,1,3)
    plot(best_GMloadingList, GMscoresAcrossGMbest, '-', 'LineWidth', LW)
    hold on
    plot(best_GMloadingList, M3scoresAcrossGMbest, '-.', 'LineWidth', LW)
    plot(best_GMloadingList, M2scoresAcrossGMbest, '--', 'LineWidth', LW)
    
    ax = gca;
    ax.YAxis.TickLabelFormat = "%.1f";
    ax.FontSize = 14;
    ylim([0 2])
    legend('GM Score', 'M3 Score', 'M2 Score', 'location', 'best')
    xlim([best_GMloadingList_lowerLimit best_GMloadingList_upperLimit])
    xlabel('Best Team GM Load Margin [-]')
    grid on
    
    % Mass Breakdown
    figure
    labels = {'Fuselage', 'Wing', 'Tail', 'Battery', 'Motor', 'Landing Gear', 'Servos', 'Wiring', 'ESC'};
    pie(best_plane.m_emptyList/best_plane.m_empty, labels);
    title("Empty Weight Breakdown")
    % ax = gca;
    % ax.FontName = "Times New Roman";
    % ax.FontSize = 18;
    % lgd = legend(labels);
    % lgd.Location = 'eastoutside';

end



%% 3-D PLOTTING

% close all;

best_GMloading = 60;
GMscores = GMscores_USC/best_GMloading;

if max(GMscores, [], 'all') > 1
    if bonuses
        GMscores(GMscores > 1) = 1 + (1 - (1./GMscores(GMscores > 1))); % Bonus for exceeding best GM multiplier
    else
        GMscores(GMscores > 1) = 1;
    end
end
scoreArray = m2scores + m3scores + GMscores;


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
    for p = 1:length(m2payloadfractionList)-4 % subtract here to not plot all payload fractions
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
            if units == "metric"
                surf(span, AR, m2payloadFraction, score)
            elseif units == "imperial"
                surf(span/0.3048, AR, m2payloadFraction, score)
            end
        end
    end
    view(3)
    view([-28 11])
    ax = gca;
    ax.FontSize = 14;
%     shading interp
    alpha 0.5
    titleString = strcat("USC GM Load Margin: ", num2str(GMMultiplierList(g)));
    title(titleString,'FontSize',14)
    if units == "metric"
        xlim([min(span, [], 'all') max(span, [], 'all')])
    elseif units == "imperial"
        xlim([min(span/0.3048, [], 'all') max(span/0.3048, [], 'all')])
    end
    xlabel("b [ft]",'FontSize',14)
    ax.XAxis.TickLabelFormat = "%.1f";
    ylim([min(AR, [], 'all') max(AR, [], 'all')])
    ylabel("AR [-]",'FontSize',14)
    ax.YAxis.TickLabelFormat = "%.f";
    % zlim([min(m2payloadFraction, [], 'all') max(m2payloadFraction, [], 'all')])
    zlabel("M2 Payload Fraction [-]",'FontSize',14)
    ax.ZAxis.TickLabelFormat = "%.2f";
    cb = colorbar; % create and label the colorbar
    cb.Label.String = 'Competition Score';
    cb.Label.FontSize = 14;
    % cb.Label.FontWeight = 'bold';
    cb.Limits = ([0 2.5]);
    grid on
    lighting phong
%         set('edgecolor',[0 0 0.4],'meshstyle','both','linewidth',.15);

end

