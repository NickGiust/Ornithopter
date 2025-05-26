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
% plane.GMloadingMultiplier = 10;
% plane.payload2Fraction = 0.35;
% plane.lengthPayload3 = 0.8;
plane.taperRatio = 0.80;

%% Define Environment
environment = defineEnvironment("Tucson");

%% 5D Trade Study
% m2payloadfractionList = 0.3:0.1:0.5;
% ARList = 4:0.5:10;
% spanList = 1.3:0.05:1.8;
% antennaLengthList = 0.6:0.1:1.2;
% GMMultiplierList = 70:10:150;

run = "M2";

% Prop Study
propStudyPlanes = [];
propDiams = 10:14;
propPitches = 6:12;

%Single plane
if run == "M2"
    %m2payloadfractionList = [0.305 0.3355];
    m2payloadfractionList = 0.3;
elseif run == "M3"
    m2payloadfractionList = 0.3;
end

ARList = 1.5/0.227;
spanList = 1.5;
antennaLengthList = 0.65;
GMMultiplierList = 70;

for i = 1:length(propDiams)
    diameter = propDiams(i);
    for j = 1:length(propPitches)
    pitch = propPitches(j);
    propellerMAT = matfile('myFile2.mat');
    allProps = propellerMAT.allProps;
    diamIndices = find(allProps(:,1)==diameter);
    pitchIndices = find(allProps(:,2)==pitch);
    [index,pos]=intersect(diamIndices,pitchIndices);
    if isempty(diamIndices) || isempty(pitchIndices) || isempty(index)
        fprintf("%i by %i does not exist\n",diameter,pitch)
        continue
    end
    diameter
    pitch
    if run == "M2"
        plane.D2 = diameter*0.0254;
        plane.P2 = pitch;
        plane.A2 = allProps(index,3);
        plane.B2 = allProps(index,4);
        plane.C2 = allProps(index,5);
        plane.E2 = allProps(index,6);
        plane.F2 = allProps(index,7);
        plane.G2 = allProps(index,8);
    end

    if run == "M3"
        plane.D3 = diameter*0.0254;
        plane.P3 = pitch;
        plane.A3 = allProps(index,3);
        plane.B3 = allProps(index,4);
        plane.C3 = allProps(index,5);
        plane.E3 = allProps(index,6);
        plane.F3 = allProps(index,7);
        plane.G3 = allProps(index,8);
    end

[m2payloadfractionArray,spanArray,antennaLengthArray,ARArray,GMMultiplierArray] = ...
    ndgrid(m2payloadfractionList,spanList,antennaLengthList,ARList, GMMultiplierList);

m2scoreArray = spanArray;
m3scoreArray = spanArray;
planes = {};
boxes = {};
m2scores_USC = spanArray;
m3scores_USC = spanArray;
GMscores_USC = spanArray;
simulationStats = [];
simulationStats.M2TOFLFails = 0;
simulationStats.M3TOFLFails = 0;
simulationStats.invalidTubes = 0;
simulationStats.boxFails = 0;

counter = 0;
totalT = 0;
totalCombos = length(m2payloadfractionList)*length(spanList)*length(antennaLengthList)*length(GMMultiplierList)*length(ARList);

tic
for p = 1:length(m2payloadfractionList)
    plane.payload2Fraction = m2payloadfractionList(p);
    
    for s = 1:length(spanList)
        plane.b = spanList(s);

        for l = 1:length(antennaLengthList)
            plane.lengthPayload3 = antennaLengthList(l);

            for a = 1:length(ARList)
                plane.AR = ARList(a);
                plane.c = plane.b/plane.AR;
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

                plane.lh = 0.9*plane.lengthPayload3; % boom is shorter than antenna to ensure box fitment
                plane.Sh = plane.Vh*plane.S*plane.c/plane.lh;
                plane.bh = sqrt(plane.ARh*plane.Sh);
                plane.ch = plane.Sh./plane.bh;
                plane.lv = plane.lh;
                plane.Sv = plane.Vv*plane.S*plane.b/plane.lv;
                plane.cv = plane.ch;
                plane.bv = plane.Sv/plane.cv;
                plane.ARv = plane.bv.^2./plane.Sv;
                % Tail
                plane.cTail = plane.ch;
                plane.bTail = plane.bh;
                plane.hTail = plane.bv;
                                
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
                    outerSectionSpanMinimumFraction = 0.83;
                    if outerSectionSpan < outerSectionSpanMinimumFraction * semispan
                        outerSectionSpan = outerSectionSpanMinimumFraction * semispan;
                    end
                    plane.outerSectionSpan = outerSectionSpan;

                    [plane, box] = checkBoxSizing(plane);

                    plane.validTubes = true;
                    [m2results, m3results, GMresults, plane] = TS_CompStrat_SimulateCompetition(plane, environment, planeToolsDirectory);
                    
                    plane.m2results = m2results;
                    plane.m3results = m3results;
                    plane.GMresults = GMresults;
                    GMMultiplierArray(p,s,l,a,g) = plane.GMloadingMultiplier;
                    
                    planes{p,s,l,a,g} = plane;
                    boxes{p,s,l,a,g} = box;
                    m2scores_USC(p,s,l,a,g) = plane.m2results.score;
                    m3scores_USC(p,s,l,a,g) = plane.m3results.score;
                    GMscores_USC(p,s,l,a,g) = plane.GMresults.score;

                    
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

                    %{
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
                    %}
                end
            end
        end
    end
end


%% Best Plane Assumptions, Print Results


%diary('14x12_6S_14x12_10S.txt');

%{
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
%}

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

best_GMloadingList_lowerLimit = 50;
best_GMloadingList_upperLimit = 50;
best_GMloadingList_increment = 5;
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
    best_plane.m2score = m2score;
    m3score = m3scores(bestScoreIndex);
    best_plane.m3score = m3score;
    GMscore = GMscores_USC(bestScoreIndex)/GMmax;

    if GMscore > 1
        if bonuses
            GMscore = 1 + (1 - (1/GMscore)); % Bonus for exceeding best GM multiplier
        else
            GMscore = 1;
        end
    end

    best_plane.GMscore = GMscore;
    
    totalscore = m2score + m3score + GMscore;

    % print results
    
    %{
    fprintf("\nFor %.d <= GMbest <= %.d:\n", GMmin, GMmax)
    fprintf("--BEST COMPETITION PLANE--\n")
    
    fprintf("\n--SCORES--\n")
    fprintf("Total Score: %.2f\n", totalscore);
    fprintf("   GM Score: %.2f\n", GMscore);
    fprintf("   M2 Score: %.2f\n", m2score); 
    fprintf("   M3 Score: %.2f\n", m3score); 
    %}
    
    fprintf("\n--SCORES--\n")
    fprintf("Total Score: %.2f\n", totalscore);
    fprintf("   GM Score: %.2f\n", GMscore);
    fprintf("   M2 Score: %.2f\n", m2score); 
    fprintf("   M3 Score: %.2f\n", m3score); 

    fprintf("\n--M2 TAKEOFF W/ %.f%% THROTTLE (MAX 60 FT)--\n", 100*best_plane.m2results.minTOthrottle)
    fprintf("M2 TOFL: %s\n", best_plane.m2results.TOFLbool)
    if best_plane.m2results.TOFLbool == 'PASSED'
        fprintf("         %.2f ft\n", best_plane.m2results.TOFL/0.3048);
        % fprintf("M2 Min. TO Throttle: %.f%%\n", 100*best_plane.m2results.minTOthrottle)
        TOsegmentIndex = find(best_plane.m2results.performance.V > best_plane.m2results.VTO, 1);
        fprintf("M2 Avg. Current, TO: %.2f A for %.2f s\n", mean(best_plane.m2results.performance.I(1:TOsegmentIndex-1)), best_plane.m2results.coursePoints.startTimes(2));
    end
    fprintf("\n--M3 TAKEOFF W/ %.f%% THROTTLE (MAX 60 FT)--\n", 100*best_plane.m3results.minTOthrottle)
    fprintf("M3 TOFL: %s\n", best_plane.m3results.TOFLbool)
    if best_plane.m3results.TOFLbool == 'PASSED'
        fprintf("         %.2f ft\n", best_plane.m3results.TOFL/0.3048);
        TOsegmentIndex = find(best_plane.m3results.performance.V > best_plane.m3results.VTO, 1);
        fprintf("M3 Avg. Current, TO: %.2f A for %.2f s\n\n", mean(best_plane.m3results.performance.I(1:TOsegmentIndex-1)), best_plane.m3results.coursePoints.startTimes(2));
    end

    if best_plane.m2results.TOFLbool == 'PASSED' && best_plane.m3results.TOFLbool == 'PASSED'
    fprintf("M2 Max Prop Torque: %.2f Nm\n", max(best_plane.m2results.performance.proptorque));
    fprintf("M3 Max Prop Torque: %.2f Nm\n", max(best_plane.m3results.performance.proptorque));
    end

    if best_plane.m2results.minTOthrottle < 1 && best_plane.m2results.TOFLbool == "PASSED"
        fprintf("\n--TAKEOFF W/ 100%% THROTTLE (MAX 60 FT)--\n")
        fprintf("M2 TOFL: %.2f ft\n", best_plane.m2TOFL_fullThrottle/0.3048);
        fprintf("M2 Avg. Current, TO: %.2f A for %.2f s\n", best_plane.m2TOFL_avgCurrent, best_plane.m2TOFL_avgCurrentTime);
    end

    if best_plane.m3results.minTOthrottle < 1 && best_plane.m3results.TOFLbool == "PASSED"
        fprintf("M3 TOFL: %.2f ft\n", best_plane.m3TOFL_fullThrottle/0.3048);
        fprintf("M3 Avg. Current, TO: %.2f A for %.2f s\n", best_plane.m3TOFL_avgCurrent, best_plane.m3TOFL_avgCurrentTime);
    end

    if best_plane.m2results.TOFLbool == "PASSED"
    fprintf("\n--LAPS/PERFORMANCE--\n")
    fprintf("M2: %.d laps in %i minutes %.1f seconds\n    %.1f%% battery energy used\n    %.f%% throttle (optimized)\n", best_plane.m2results.nLaps, floor(best_plane.m2results.missionTime/60),rem(best_plane.m2results.missionTime,60),100*best_plane.m2results.batteryUsage,100*best_plane.m2results.optimizedThrottleSetting);
    end
    if best_plane.m3results.TOFLbool == "PASSED"
    fprintf("M3: 3 laps in %.1f seconds\n    %.1f%% battery energy used\n", best_plane.m3results.missionTime, 100-100*best_plane.m3results.performance.s(end));
    end

    if best_plane.m2results.TOFLbool == 'PASSED' && best_plane.m3results.TOFLbool == 'PASSED'
    fprintf("M2 Max. Cruise Speed: %.0f mph\n", max(best_plane.m2results.performance.V)*2.237)
    fprintf("M3 Max. Cruise Speed: %.0f mph\n", max(best_plane.m3results.performance.V)*2.237)
    end

    %{
    fprintf("\n--PAYLOAD--\n")
    fprintf("M2 Payload Fraction: %.2f\n", best_plane.payload2Fraction);
    fprintf("M2 Payload Mass: %.2f kg\n", best_plane.payload2Fraction*best_plane.m2);
    fprintf("M3 Antenna Length: %.2f m\n", best_plane.lengthPayload3);
    fprintf("GM Loading Multiplier (Design Point): %.2f\n", best_plane.GMloadingMultiplier)
    fprintf("GM Loading Multiplier (True): %.2f\n", best_plane.trueGMloadingMultiplier)
    fprintf("GM Load: %.f kg\n", best_plane.trueGMloadingMultiplier*best_plane.m2)
    %}
    fprintf("\n--MASS--\n")
    fprintf("Empty Mass: %.2f kg\n", best_plane.m_empty);
    fprintf("M2 Mass: %.2f kg\n", best_plane.m2);
    fprintf("M2 Payload Fraction: %.2f\n", best_plane.mPayload2/best_plane.m2);
    fprintf("M2 Payload Mass: %.2f kg\n", best_plane.mPayload2);
    fprintf("M3 Mass: %.2f kg\n", best_plane.m3);
    fprintf("Wing Mass: %.2f kg\n", best_plane.mWing);
    %{
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
    %}
    
end

    propStudyPlanes = [propStudyPlanes best_plane];
    end
end

%% Excel File

propStudyData = [];
if run == "M2"
    for i = 1:length(propStudyPlanes)
        best_plane = propStudyPlanes(i);
        if best_plane.m2results.TOFLbool == 'PASSED' && best_plane.m2results.nLaps < 30
            propStudyData(i,1) = round(best_plane.D2 * 39.37);
            propStudyData(i,2) = best_plane.P2;
            propStudyData(i,3) = best_plane.m2score;
    
            propStudyData(i,4) = best_plane.m2results.minTOthrottle;
            propStudyData(i,5) = best_plane.m2results.TOFL/0.3048;
            TOsegmentIndex = find(best_plane.m2results.performance.V > best_plane.m2results.VTO, 1);
            propStudyData(i,6) = mean(best_plane.m2results.performance.I(1:TOsegmentIndex-1));
            propStudyData(i,7) = max(best_plane.m2results.performance.proptorque);

            if best_plane.m2results.minTOthrottle < 1 && best_plane.m2results.TOFLbool == "PASSED"
                propStudyData(i,8) = 1;
                propStudyData(i,9) = best_plane.m2TOFL_fullThrottle/0.3048;
                propStudyData(i,10) = best_plane.m2TOFL_avgCurrent;
            else
                propStudyData(i,8) = 1;
                propStudyData(i,9) = best_plane.m2results.TOFL/0.3048;
                propStudyData(i,10) = mean(best_plane.m2results.performance.I(1:TOsegmentIndex-1));
            end

            propStudyData(i,11) = best_plane.m2results.nLaps;
            propStudyData(i,12) = best_plane.m2results.batteryUsage;
            propStudyData(i,13) = best_plane.m2results.optimizedThrottleSetting;
            propStudyData(i,14) = max(best_plane.m2results.performance.V)*2.237;
            propStudyData(i,15) = max(best_plane.m2results.Vstall)*2.237;

            propStudyData(i,16) = best_plane.payload2Fraction;
            propStudyData(i,17) = best_plane.payload2Fraction*best_plane.m2;
        end
    end
end
writematrix(propStudyData,'M2compStrat23.xlsx', 'Sheet', 'M2 - 3MPS', 'Range', 'B2')

%% PROP STUDIES

% Plot prop study results
figure
legentries = [];
markerStyles = {'o', 's', 'd', 'p', 'h', 'v', '^', '<', '>', 'x'};
if run == "M2"
    for i = 1:length(propStudyPlanes)
        planefrompropstudy = propStudyPlanes(i);
        if planefrompropstudy.m2results.TOFLbool == 'PASSED' && planefrompropstudy.m2results.nLaps < 30
            marker = markerStyles{mod(i-1, numel(markerStyles)) + 1};
            s = scatter3(planefrompropstudy.m2results.TOFL*3.281,planefrompropstudy.mPayload2,planefrompropstudy.m2score+planefrompropstudy.GMscore,250, marker, 'filled');
            legentries = [legentries strcat(string(round(planefrompropstudy.D2*39.3701)),"x",string(planefrompropstudy.P2))];
        end
        hold on
    end
    title(strcat(run," Prop Study w/ ",string(planefrompropstudy.nSeries2),"S (",string(environment.windSpeed), " m/s Wind)"))
    xlabel("TOFL [ft]",'FontSize',18,'FontWeight','bold')
    ylabel("Best M2 Payload [kg]",'FontSize',18,'FontWeight','bold')
    zlabel("M2 + GM Score [-]",'FontSize',18,'FontWeight','bold')
    %ylim([0.3 0.5])
    zlim([1 2])
    ax = gca;
    ax.FontSize = 14;
    ax.ZTick = unique( round(ax.ZTick) );
    l = legend(legentries);
    l.Location = 'bestoutside';
    l.FontSize = 24;
    view([10 10])
    grid on
end

if run == "M3"
    for i = 1:length(propStudyPlanes)
        planefrompropstudy = propStudyPlanes(i);
        if planefrompropstudy.m3results.TOFLbool == 'PASSED'
            marker = markerStyles{mod(i-1, numel(markerStyles)) + 1};
            s = scatter3(planefrompropstudy.m3results.TOFL*3.281,max(planefrompropstudy.m3results.performance.proptorque),planefrompropstudy.m3results.missionTime,250, marker, 'filled');
            legentries = [legentries strcat(string(round(planefrompropstudy.D3*39.3701)),"x",string(planefrompropstudy.P3))];
        end
        hold on
    end
    
    title(strcat(run," Prop Study w/ ",string(planefrompropstudy.nSeries3),"S (",string(environment.windSpeed), " m/s Wind)"))
    xlabel("TOFL [ft]",'FontSize',18,'FontWeight','bold')
    ylabel("Torque [Nm]",'FontSize',18,'FontWeight','bold')
    zlabel("Time [s]",'FontSize',18,'FontWeight','bold')
    ax = gca;
    ax.FontSize = 14;
    l = legend(legentries);
    l.Location = 'north';
    l.FontSize = 18;
    view([10 10])
    grid on
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
