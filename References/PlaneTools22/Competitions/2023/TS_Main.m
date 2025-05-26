clc; clear; close all;
global planeToolsDirectory;

%% Inputs
currentDirectory = split(pwd, 'PlaneTools22');
planeToolsDirectory = [currentDirectory{1} '\PlaneTools22'];

%% Read aircraft from text file
plane = definePlane()
plane.GMloadingMultiplier = 10;
plane.payload2Fraction = 0.35;
plane.lengthPayload3 = 0.8;
plane.taperRatio = 0.7;
%% 3D Trade Study

% SPAN VS AR LIST VS GM LOADING
spanList = 0.8:0.1:1;
m2payloadfractionList = 0.3:0.05:0.4;
M2lapsList = 0; 
[X,Y,Z] = meshgrid(spanList,m2payloadfractionList,M2lapsList);
scoreList = X;
m2scoreList = X;
m3scoreList = X;
planes = {};
boxes = {};
counter = 0;
best_finalscore = 0;
totalCombos = length(spanList)*length(m2payloadfractionList)*length(M2lapsList);
for i = 1:length(spanList)
    for j = 1:length(m2payloadfractionList)
            plane.b = spanList(i);
            plane.payload2Fraction = m2payloadfractionList(j);
            plane.c = plane.b/plane.AR;
            [planeFits_boolean, box, plane] = checkBoxSizing(plane);
            if planeFits_boolean == true
                [m2results, m3results, finalscore, plane] = TS_SimulateCompetition(plane, planeToolsDirectory);
                Z(j,i) = m2results.nLaps;
                plane.m2results = m2results;
                plane.m3results = m3results;
                scoreList(j,i) = finalscore;
                m2scoreList(j,i) = m2results.score;
                m3scoreList(j,i) = m3results.score;
            else
                scoreList(j,i) = 0;
                m2scoreList(j,i) = 0;
                m3scoreList(j,i) = 0;
            end
            counter = counter + 1;
            %clc
            fprintf("Total # of Planes in Trade Study: %i\n", totalCombos)
            fprintf("Simulation %.2f%% Completed...\n", (counter/totalCombos)*100)
            planes{j,i} = plane;
            boxes{j,i} = box;
    end
end

% % SPAN VS M2 PAYLOAD FRACTION VS GM LOADING
% plane.AR = 6;
% spanList = 1:0.1:2;
% antennaLengthList = 0.5:0.02:1;
% M3timeList = linspace(10,20,1);
% [X,Y,Z] = meshgrid(spanList,antennaLengthList,M3timeList);
% scoreList = X;
% m2scoreList = X;
% m3scoreList = X;
% planes = {};
% boxes = {};
% counter = 0;
% best_finalscore = 0;
% totalCombos = length(spanList)*length(antennaLengthList)*length(M3timeList);
% for i = 1:length(spanList)
%     for j = 1:length(antennaLengthList)
%         for k = 1:length(M3timeList)
%             plane.b = spanList(i);
%             plane.hWinglet = plane.b*0.15;
%             plane.lengthPayload3 = antennaLengthList(j);
%             plane.c = plane.b/plane.AR;
%             [planeFits_boolean, box] = checkBoxSizing(plane);
%             if planeFits_boolean == true
%                 [m2results, m3results, finalscore, plane] = TS_SimulateCompetition(plane, planeToolsDirectory);
%                 Z(j,i,k) = m3results.missionTime;
%                 plane.m2results = m2results;
%                 plane.m3results = m3results;
%                 scoreList(j,i,k) = finalscore;
%                 m2scoreList(j,i,k) = m2results.score;
%                 m3scoreList(j,i,k) = m3results.score;
%             else
%                 scoreList(j,i,k) = 0;
%                 m2scoreList(j,i,k) = 0;
%                 m3scoreList(j,i,k) = 0;
%                 Z(j,i,k) = 300;
%             end
%             counter = counter + 1;
%             %clc
%             fprintf("Total # of Planes in Trade Study: %i\n", totalCombos)
%             fprintf("Simulation %.2f%% Completed...\n", (counter/totalCombos)*100)
%             planes{j,i,k} = plane;
%             boxes{j,i,k} = box;
%         end
%     end
% end


% antennaLengthList = 0.5:0.1:1;
% m2PayloadFractionList = 0.35:0.1:0.75;
% GMloadingList = 10:10:100; % multiplier on M2 mass
% [X,Y,Z] = meshgrid(antennaLengthList,m2PayloadFractionList,GMloadingList);
% scoreList = X;
% m2scoreList = X;
% m3scoreList = X;
% planes = {};
% boxes = {};
% counter = 0;
% best_finalscore = 0;
% totalCombos = length(antennaLengthList)*length(m2PayloadFractionList)*length(GMloadingList);
% for i = 1:length(antennaLengthList)
%     for j = 1:length(m2PayloadFractionList)
%         for k = 1:length(GMloadingList)
%             plane.lengthPayload3 = antennaLengthList(i);
%             plane.DApayload3 = (1.17*(plane.lengthPayload3*0.0214))/plane.S;
%             plane.hWinglet = plane.b*0.15;
%             plane.payload2Fraction = m2PayloadFractionList(j);
%             plane.GMloadingMultiplier = GMloadingList(k);
%             [planeFits_boolean, box] = checkBoxSizing(plane);
%             if planeFits_boolean == true
%                 [m2results, m3results, finalscore, plane] = TS_SimulateCompetition(plane, planeToolsDirectory);
%                 Z(j,i,k) = plane.GMloadingMultiplier;
%                 plane.m2results = m2results;
%                 plane.m3results = m3results;
%                 scoreList(j,i,k) = finalscore;
%                 m2scoreList(j,i,k) = m2results.score;
%                 m3scoreList(j,i,k) = m3results.score;
%             else
%                 scoreList(j,i,k) = 0;
%                 m2scoreList(j,i,k) = 0;
%                 m3scoreList(j,i,k) = 0;
%             end
%             counter = counter + 1;
%             %clc
%             fprintf("Total # of Planes in Trade Study: %i\n", totalCombos)
%             fprintf("Simulation %.2f%% Completed...\n", (counter/totalCombos)*100)
%             planes{j,i,k} = plane;
%             boxes{j,i,k} = box;
%         end
%     end
% end

%% 3D plotting
close all;

% % 3D trade
% figure;
% h = scatter3(X(:),Y(:),Z(:), [], scoreList(:),'filled'); 
% colorbar
% view(-80,15)
% set(h, 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', 0.5)
% grid on   
% %xlabel("M3 Antenna Length [m]")
% xlabel("Span [m]", 'FontSize',12,'FontWeight','bold')
% %ylabel("M2 Payload Fraction")
% ylabel("AR",'FontSize',12,'FontWeight','bold')
% %zlabel("GM Loading Multiplier")
% zlabel("GM Load Margin [kg]",'FontSize',12,'FontWeight','bold')

figure
hold on
for i = 1:size(X,3)
    s1 = surf(X(:,:,i),Y(:,:,i),Z(:,:,i),m2scoreList(:,:,i))
end
view(3)
shading interp
alpha 0.5
xlabel("Span [m]",'FontSize',12,'FontWeight','bold')
ylabel("M2 Payload Fraction [-]",'FontSize',12,'FontWeight','bold')
zlabel("M2 Laps",'FontSize',12,'FontWeight','bold')
cb = colorbar;                                  % create and label the colorbar
cb.Label.String = 'M2 Score';
cb.Label.FontSize=12;
grid on
lighting phong
set(s1,'edgecolor',[0 0 0.4],'meshstyle','both','linewidth',.15);

% BEST OVERALL SCORING PLANE
[val,idx] = max(scoreList(:)) ;  % get maximum in Z
best_plane = planes(idx);
best_plane = best_plane{1};
best_box = boxes(idx);
best_box = best_box{1};

fprintf("\n\n\n--BEST COMPETITION PLANE--\n")

fprintf("\n--SCORES--\n")
fprintf("Total Score: %.2f\n", val);
% fprintf("GM Score: %.2f\n", best_plane.GMloadingMultiplier/100);
fprintf("GM Score: %.2f\n", val-best_plane.m2results.score-best_plane.m3results.score);
fprintf("M2 Score: %.2f\n", best_plane.m2results.score); 
fprintf("M3 Score: %.2f\n", best_plane.m3results.score); 

fprintf("\n--TOFL (MAX 60 FT)--\n")
fprintf("M2 TOFL: %s, %.2fft\n", best_plane.m2results.TOFLbool, best_plane.m2results.TOFL*3.281); 
fprintf("M3 TOFL: %s, %.2fft\n", best_plane.m3results.TOFLbool, best_plane.m3results.TOFL*3.281);

fprintf("\n--LAPS--\n")
fprintf("M2 LAPS: %.f\n", best_plane.m2results.nLaps); 
fprintf("M3 TIME (3 LAPS): %.2f sec\n", best_plane.m3results.missionTime);

fprintf("\n--PAYLOAD--\n")
fprintf("M2 Payload Fraction: %.2f\n", best_plane.payload2Fraction);
fprintf("M3 Antenna Length: %.2fm\n", best_plane.lengthPayload3);

fprintf("\n--MASS--\n")
fprintf("M2 Mass: %.2fkg\n", best_plane.m2);
fprintf("M3 Mass: %.2fkg\n", best_plane.m3);
fprintf("Wing Mass: %.2fkg\n", best_plane.mWing);
fprintf("Wing # of Spar Cap Plies (one side): %.f\n", max(best_plane.wingSparCapPlies));
fprintf("GM Loading Multiplier: %.2f\n", best_plane.GMloadingMultiplier)
fprintf("GM Load: %.fkg\n", best_plane.GMloadingMultiplier*best_plane.m2)

fprintf("\n--WING--\n")
fprintf("Span: %.2fm\n", best_plane.b)
fprintf("Chord: %.2fm\n", best_plane.c)
fprintf("AR: %.f\n", best_plane.AR)

fprintf("\n--TAIL--\n")
fprintf("Span: %.2fm\n", best_plane.bTail)
fprintf("Chord: %.2fm\n", best_plane.cTail)
fprintf("Height: %.2fm\n", best_plane.hTail)

fprintf("\n--BOX--\n")
fprintf("Box length: %.2fin\n", best_box.length*39.37)
fprintf("Box width: %.2fin\n", best_box.width*39.37)
fprintf("Box height: %.2fin\n", best_box.height*39.37)
fprintf("Sum of Dimensions = %.2fin\n", (best_box.length+best_box.width+best_box.height)*39.37);

% BEST M2 SCORING PLANE
[val,idx] = max(m2scoreList(:)) ;  % get maximum in Z 
fprintf("\n--Best M2 Plane--\n")
fprintf("Span: %.2fm\n", X(idx))
fprintf("AR: %.f\n", Y(idx))

% BEST M3 SCORING PLANE
[val,idx] = max(m3scoreList(:)) ;  % get maximum in Z 
fprintf("\n--Best M3 Plane--\n")
fprintf("Span: %.2fm\n", X(idx))
fprintf("AR: %.f\n", Y(idx))
%% Results
% [m2results, m2mass, m3results, m3mass, finalscore] = TS_SimulateCompetition(plane, planeToolsDirectory);
% printMissionResults(m2results, m2mass);
% printMissionResults(m3results, m3mass);
% fprintf('---------------- FINAL SCORE: %.2f ----------------\n', finalscore)
