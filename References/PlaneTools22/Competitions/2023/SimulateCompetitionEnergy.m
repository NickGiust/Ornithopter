clc;
clear;
warning('off','all');

%% PLANETOOLS 2022 ~~ RUBBER 2.0 ~~
global planeToolsDirectory;
year = '2022';
seedPlane = 'CDRplane';
location = 'Wichita';

%% Add directories to path and initialize global propulsion functions
addpath([planeToolsDirectory '\PhysicsModels\']);
addpath([planeToolsDirectory '\ComponentLibrary']);
addpath([planeToolsDirectory '\Plotting']);
addpath([planeToolsDirectory '\Environments']);

%% Aircraft
pl = ReadPlane(seedPlane, year, planeToolsDirectory);

%% Environment parameters
dens = 1.225;

%% Estimated parameters ("fudge" factors)
% Proportion of static thrust produced at max efficiency
dynamicShare = .35;
% Value during turns
dynamicShareTurn = .5;
% Total propulsive efficiency
eta = .65;
% Stall thrust fractions of theoretical static thrust
TstallFrac2 = 1;
TstallFrac3 = 1;
% Maximum CL
CLmax = 2.5;
% CD for spoilers
CDspoiler = 1.1;
% Headwind
headwind = 5;
% Braking resistance on stopping
pl.br = 0.8;
% Taxi distance
xTaxi = 5;
% Taxi speed
Vtaxi = 2;
% Power consumption during deceleration, stopping, taxi, deployment
PelLow = .15;
% Electrical power limit
Plimit = 20000;
% Deployment time
tDeployPackage = 2.5;
% Ground mission times
tLoadSyringe = .3;
tUnloadSyringe = .1;
tLoadPackage = 5;
tRun = 12;
% Technological advancement (total ability to reduce weight and drag)
tech = 1;

TOFLlimit = 7.62;

%% Scoring Assumptions
syringeSpeedBest = 10;
nLapsBest = 15;
tGMbest = 10;

%% TRADE STUDY PARAMETERS

%chords = .4;
chords = .2:.3:.8;
%chords = .2:.04:1.2;

%noPackages = 9;
noPackages = 1:2:12;
%noPackages = 1:1:14;

%noSyringes = 470;
noSyringes = [50 100 400 600];
%noSyringes = 10:20:600;

%staticThrust3 = 70;
staticThrust3 = [20 40 60 100 120 150];
%staticThrust3 = 40:10:80;

%staticThrust2 = 158;
staticThrust2 = [150 200 250 300 350];
%staticThrust2 = 200:20:300;

%spans = 2.44;
spans = [1.4 1.8 2 2.44];
%spans = 2.44;


bestScores = [];
bestSpans = [];
bestChords = [];
bestNP2 = [];
bestNP3 = [];
bestTS2 = [];
bestTS3 = [];
    
bestScore = 0;
planeNo = 0;
noScoringPlanes = 0;
scores = zeros(3872,4);

Wtech = 1;
Dtech = 1;
eTech = 1;

percent = 1;
normal = 4;
%varyRange = [1.5 1.8 2 2.5 3 3.5];
varyRange = 2.3;
category = 'Maximum lift coefficient (TO & turns)';
variable = 'C_{L,max}';

%% Score sensitivity parameter
for scoreParam = varyRange
    CLmax = scoreParam;
    
    bestScore = 0;
    %planeNo = 0;
    
    %% Chord
    for wingspan = spans
        pl.b = wingspan;

        for chord = chords
            pl.c = chord;

            %% Payloads
            for nPayloads3 = noPackages
                pl.nPayloads3 = nPayloads3;
                
                for nPayloads2 = noSyringes
                    pl.nPayloads2 = nPayloads2;
                    
                    if pl.nPayloads2 < 10*nPayloads3
                        continue;
                    end

                    for Tstatic2 = staticThrust2

                        for Tstatic3 = staticThrust3
       
                            planeNo = planeNo + 1
                            
                            %% Fuselage, Landing Gear, and Tail Adjustments
                            S = pl.b*pl.c;
                            Cvt = .08;
                            Cht = .08;
                            pl.lFuse = max([(pl.b/2) (.3+(.02*pl.nPayloads2))]);
                            Svt = Cvt*pl.b/(pl.lFuse*.7);
                            Sht = Cht*pl.b/(pl.lFuse*.7);
                            pl.cTail = (Sht/2.5)^.5;
                            pl.bTail = 2.5*pl.cTail;
                            pl.hTail = Svt/pl.cTail;
                            pl.hWing = .2*pl.b;

                            if pl.lFuse > 2.44
                                pl.lFuse = 2.44;
                                pl.hFuse = .05 + (.03 * (pl.nPayloads2/(11*5)));
                            end

                            %% Mass and drag buildups
                            m = massBuildup(pl, 2).m*tech*Wtech;
                            CD0 = dragBuildup(pl, 2).CD0*tech*Dtech;
                            e = dragBuildup(pl, 3).e*eTech;
                            W = m*9.81;

                            AR = (pl.b^2)/S;
                            k = 1/(pi*e*AR);
                            S = (pl.b-pl.wFuse)*pl.c;

                            %% MISSION 1
                            score1 = 1;

                            %% MISSION 2
                            %% Cruise (speed)
                            Vcruise = (((dynamicShare*Tstatic2)+((((dynamicShare^2)*(Tstatic2^2))-(4*k*(W^2)*CD0))^.5))/(CD0*dens*S))^.5;
                            Pel = (dynamicShare*Tstatic2*Vcruise)/eta;
                            Vcruise2 = Vcruise;
                            Pel2 = Pel;

                            if ~isreal(Vcruise)
                                continue;
                            end

                            %% Turn (time)
                            % Increased wing area from leading edge slats
                            S = S*1.2;
                            %Vturn = ((2*dynamicShareTurn*Tstatic2)/(dens*S*(CD0+(k*(CLmax^2)))))^.5;
                            Vturn = .7*Vcruise;
                            Lturn = 0.5*dens*(Vturn^2)*S*CLmax;
                            loadFactor = Lturn/W;
                            S = S/1.2;
                            
                            if loadFactor > pl.nStruct
                                loadFactor = pl.nStruct;
                            end
                            
                            r = (Vturn^2)/(9.81*loadFactor);
                            d180 = pi*r;
                            t180 = d180/Vturn;
                            
                            loadFactor2 = loadFactor;

                            %% Takeoff (distance & time)
                            m2 = m;
                            Vstall = sqrt((2*W)/(CLmax*dens*S));
                            VTO = 1.2*Vstall;
                            gfx = ((16*pl.hWing/pl.b)^2)/(1+((16*pl.hWing/pl.b)^2));

                            % Acceleration = force/mass
                            Tstall2 = TstallFrac2*Tstatic2;
                            staticThr = Tstall2/m;
                            drag = (dens*S*(CD0 + (k*gfx*(pl.CLgroundRoll^2)) - (pl.rr*pl.CLgroundRoll)))/(2*m);
                            fric = (pl.rr*W)/m;

                            % V as a function of t up to takeoff speed
                            tspan = 0:.1:20;
                            [t,V] = ode45((@(t,V) (staticThr*(1-((1-dynamicShare)*(V/Vcruise)))) ...
                                -fric-(drag*(V^2))), tspan, headwind);

                            TOindex = find( V > VTO, 1 );

                            % Distance as a function of t up to takeoff speed
                            try
                                x = cumsum(diff(t(1:(TOindex+1))).*(V(1:(TOindex))-headwind));
                            catch
                                x = 100;
                            end

                            % Takeoff field length
                            if isempty(x) || x(end) > TOFLlimit
                                continue;
                            else
                                TOFL2 = x(end);
                                tTO = t(TOindex);
                            end

                            %% Score
                            % Time to fly a lap
                            tLap = (612/Vcruise) + (t180*4);

                            % Add wind drift
                            tWindDrift = (tLap*headwind)/Vcruise;
                            tLap = tLap + tWindDrift;

                            tM2 = 3*tLap + (2*tTO);

                            if Pel > Plimit
                                continue;
                            end

                            % Energy to fly a lap
                            Elap = Pel*tLap;

                            % No of laps to 100 Wh
                            nLapsEnergy = floor(100*3600/Elap);

                            % Score 0 if less than 3 laps are completed
                            if nLapsEnergy < 3
                                score2 = 0;
                            else
                                score2 = 1+((pl.nPayloads2/tM2)/syringeSpeedBest);
                            end

                            %% MISSION 3
                            m = massBuildup(pl, 3).m*tech*Wtech;
                            W = m*9.81;

                            %% Cruise (speed)
                            Wcruise = W - (.5*pl.mPayload3*pl.nPayloads3);
                            Vcruise = (((dynamicShare*Tstatic3)+((((dynamicShare^2)*(Tstatic3^2))-(4*k*(Wcruise^2)*CD0))^.5))/(CD0*dens*S))^.5;
                            Pel = (dynamicShare*Tstatic3*Vcruise)/eta;
                            Vcruise3 = Vcruise;
                            Pel3  = Pel;

                            if ~isreal(Vcruise)
                                continue;
                            end

                            %% Takeoff (distance/time)
                            Vstall = sqrt((2*W)/(CLmax*dens*S));
                            VTO = 1.2*Vstall;
                            gfx = ((16*pl.hWing/pl.b)^2)/(1+((16*pl.hWing/pl.b)^2));

                            % Acceleration = force/mass
                            Tstall3 = TstallFrac3*Tstatic3;
                            staticThr = Tstatic3/m;
                            drag = (dens*S*(CD0 + (k*gfx*(pl.CLgroundRoll^2)) - (pl.rr*pl.CLgroundRoll)))/(2*m);
                            fric = (pl.rr*W)/m;

                            % V as a function of t up to takeoff speed
                            tspan = 0:.1:20;
                            [t,V] = ode45((@(t,V) (staticThr*(1-((1-dynamicShare)*(V/Vcruise)))) ...
                                -fric-(drag*(V^2))), tspan, headwind);

                            TOindex = find( V > VTO, 1 );
                            tTO = t(TOindex);

                            % Distance as a function of t up to takeoff speed
                            try
                                x = cumsum(diff(t(1:(TOindex+1))).*(V(1:(TOindex))-headwind));
                            catch
                                x = 100;
                            end
                            
                            % Takeoff field length
                            if isempty(x) || x(end) > TOFLlimit
                                continue;
                            else
                                TOFL3 = x(end);
                            end

                            %% Accelerate (distance/time)
                            % Acceleration = force/mass
                            dragParasite = (dens*S*CD0)/(2*m);
                            dragInduced = ((2*k*(W^2))/(dens*S))/m;

                            % V as a function of t up to near-cruise speed
                            tspan = 0:.1:30;
                            [t,V] = ode45((@(t,V) (staticThr*(1-((1-dynamicShare)*(V/Vcruise)))) ...
                                -(dragParasite*(V^2))-(dragInduced/(V^2))), tspan, VTO);

                            crIndex = find( V > 0.9*Vcruise, 1 );
                            tAcc = t(crIndex);

                            % Distance as a function of t up to near-cruise speed
                            x = cumsum(diff(t(1:(crIndex+1))).*(V(1:(crIndex))));

                            % Acceleration distance
                            if isempty(x)
                                continue;
                            else
                                xAcc = x(end);
                            end
                            
                            %% Turn (time)
                            % Increased wing area from leading edge slats
                            S = S*1.2;
                            %Vturn = ((2*dynamicShareTurn*Tstatic3)/(dens*S*(CD0+(k*(CLmax^2)))))^.5;
                            Vturn = .7*Vcruise;
                            Lturn = 0.5*dens*(Vturn^2)*S*CLmax;
                            loadFactor = Lturn/W;
                            S = S/1.2;
                            
                            if loadFactor > pl.nStruct
                                loadFactor = pl.nStruct;
                            end
                            
                            r = (Vturn^2)/(9.81*loadFactor);
                            d180 = pi*r;
                            t180 = d180/Vturn;

                            %% Decelerate (distance/time)
                            % Spoiler drag area
                            DAspoiler = CDspoiler*(pl.c/4)*(pl.b/2);

                            % Acceleration = force/mass
                            dragParasite = (dens*S*CD0)/(2*m);
                            dragSpoiler = (dens*DAspoiler)/(2*m);
                            dragInduced = ((2*k*(W^2))/(dens*S))/m;

                            % V as a function of t up to landing (takeoff) speed
                            tspan = 0:.1:30;
                            [t,V] = ode45((@(t,V) -((dragParasite+dragSpoiler)*(V^2))-(dragInduced/(V^2))), tspan, Vcruise);

                            landIndex = find( V < VTO, 1 );
                            tDec = t(landIndex);

                            % Distance as a function of t up to landing (takeoff) speed
                            x = cumsum(diff(t(1:(landIndex+1))).*(V(1:(landIndex))));

                            % Deceleration distance
                            xDec = x(end);

                            %% Stop (distance/time)
                            tspan = 0:.1:10;
                            drag = (dens*S*(CD0 + (k*gfx*(pl.CLgroundRoll^2)) - (pl.rr*pl.CLgroundRoll)))/(2*m);
                            fric = (pl.br*W)/m;
                            [t,V] = ode45((@(t,V) -fric-(drag*(V^2))), tspan, VTO);

                            stopIndex = find( V < 0, 1 );
                            tStop = t(stopIndex);

                            % Distance as a function of t to stop
                            x = cumsum(diff(t(1:(stopIndex+1))).*(V(1:(stopIndex))-headwind));

                            % Takeoff field length
                            xStop = x(end);

                            %% Taxi (time)
                            tTaxi = xTaxi/Vtaxi;

                            %% Scoring
                            % Add takeoff distance, acceleration distance,
                            % deceleration distance, stop distance, and taxi
                            % distance
                            xNonCruiseStraight = TOFL3 + xAcc + xDec + xStop + xTaxi;

                            % If greater than 306 m (1000 ft), this distance is the
                            % new straightaway distance
                            xStraightaway = max([306 xNonCruiseStraight]);

                            % If less, find the additional straightaway distance 
                            % needed to complete upwind leg
                            xAdditionalUpwind = 306 - xNonCruiseStraight;

                            % Find total lap time
                            tLap = tTO + tAcc + (xStraightaway/Vcruise) + (t180*4) ...
                                + (xAdditionalUpwind/Vcruise) + tDec + tStop + tTaxi ...
                                + tDeployPackage;
                            %{
                            figure();
                            pie([tTO tAcc (xStraightaway/Vcruise) (t180*4) (xAdditionalUpwind/Vcruise) ...
                                tDec tStop tTaxi tDeployPackage], ...
                                [0 1 0 0 0 0 0 0 1], ...
                                ["Takeoff", "Acceleration", "Downwind Straightaway", "Turns", ...
                                "Upwind Straightaway", "Deceleration", "Stopping", "Taxiing", "Deployment"]);
                            %}

                            % Add wind drift time
                            tWindDrift = (tLap*headwind)/Vcruise;
                            tLap = tLap + tWindDrift;

                            % No of laps to 10 min
                            nLapsTime = floor(600/tLap);

                            % Low-power share of lap time
                            tLow = tDec + tStop + tTaxi + tDeployPackage;

                            % Energy to fly a lap
                            Elap = Pel*((tLap-tLow)+(PelLow*tLow));

                            if Pel > Plimit
                                continue;
                            end

                            % No of laps to 100 Wh
                            nLapsEnergy = floor(100*3600/Elap);

                            % No of laps during which packages can be deployed
                            nLapsPayloads = pl.nPayloads3;

                            % Whichever is less is the number of laps
                            nLaps = min([nLapsTime nLapsEnergy nLapsPayloads]);

                            if nLaps < 0
                                continue;
                            elseif nLaps < nLapsTime
                                missionTime = tLap*nLaps;
                            else
                                missionTime = 600;
                            end

                            % Score mission
                            score3 = 2 + (nLaps/nLapsBest);

                            %% GROUND MISSION
                            tGM = (tLoadSyringe + tUnloadSyringe)*pl.nPayloads2 + tLoadPackage*pl.nPayloads3 + tRun;
                            if tGM > 300
                                continue;
                            end
                            scoreGM = tGMbest/tGM;

                            %% Score plane and compare
                            score = score1 + score2 + score3 + scoreGM;

                            % Compare to previous best scoring aircraft
                            % If better, save score, no of laps, no of eggs, time, chord, no of
                            % payloads, pitch angle, Kv, voltage, and diameter
                            if score > bestScore
                                bestPlane.no = planeNo;
                                bestPlane.score = score;
                                bestPlane.score2 = score2;
                                bestPlane.score3 = score3;
                                bestPlane.scoreGM = scoreGM;
                                
                                bestPlane.wingspan = pl.b;
                                bestPlane.chord = pl.c;
                                bestPlane.fuseHeight = pl.hFuse;
                                bestPlane.CD0 = CD0;
                                bestPlane.nSyringes = pl.nPayloads2;
                                bestPlane.nPackages = pl.nPayloads3;
                                bestPlane.Tstatic2 = Tstatic2;
                                bestPlane.Tstatic3 = Tstatic3;
                                
                                bestPlane.syringeSpeed = pl.nPayloads2/tM2;
                                bestPlane.nLaps3 = nLaps;
                                bestPlane.tGM = tGM;

                                bestPlane.tM2 = tM2;
                                bestPlane.TOFL2 = TOFL2;
                                bestPlane.m2 = m2;
                                bestPlane.PF2 = (pl.mPayload2*pl.nPayloads2)/m2;
                                bestPlane.Vcruise2 = Vcruise2;
                                bestPlane.CL2 = (m2*2*9.81)/(dens*(Vcruise2^2)*S);
                                bestPlane.loadFactor2 = loadFactor2;
                                bestPlane.Pel2 = Pel2;

                                bestPlane.tM3 = missionTime;
                                bestPlane.TOFL3 = TOFL3;
                                bestPlane.mTO3 = m;
                                bestPlane.PF3 = (pl.mPayload3*pl.nPayloads3)/m;
                                bestPlane.Vcruise3 = Vcruise3;
                                bestPlane.CL3 = (m*2*9.81)/(dens*(Vcruise3^2)*S);
                                bestPlane.loadFactor3 = loadFactor;
                                bestPlane.Pel3 = Pel3;

                                bestScore = score;
                            end                   
                            noScoringPlanes = noScoringPlanes + 1;
                            scores(noScoringPlanes,:) = [pl.c pl.nPayloads2 pl.nPayloads3 score];
                        end
                    end
                end
            end
        end
    end
    try
        bestSpans = [bestSpans bestPlane.wingspan];
        bestChords = [bestChords bestPlane.chord];
        bestNP2 = [bestNP2 bestPlane.nSyringes];
        bestNP3 = [bestNP3 bestPlane.nPackages];
        bestTS2 = [bestTS2 bestPlane.Tstatic2];
        bestTS3 = [bestTS3 bestPlane.Tstatic3];
    catch
        bestSpans = [bestSpans 0];
        bestChords = [bestChords 0];
        bestNP2 = [bestNP2 0];
        bestNP3 = [bestNP3 0];
        bestTS2 = [bestTS2 0];
        bestTS3 = [bestTS3 0];
    end
    bestScores = [bestScores bestScore];
end

noPlanes = planeNo
noScoringPlanes
bestPlane

heat = scores(:,4)-min(scores(:,4));
heat = heat/max(heat);

heatR = 2*(heat-.5);
heatR(heatR < 0) = 0;

heatB = 2*(.5-heat);
heatB(heatB < 0) = 0;

heatG = 1-heatR-heatB;

scatter3(scores(:,1), scores(:,2), scores(:,3),[], [heatR heatG heatB]);

%{
hold on;

plot(varyRange*percent, bestScores, 'LineWidth', 4.5);
xlabel(variable);
title(category);
legend('Score', 'Span', 'Chord', 'N_{payloads,2}', ...
'N_{payloads,3}', 'T_{static,2}', 'T_{static,3}');
ylabel('Best DBF Score');
ylim([3 7]);


yyaxis right;

plot(varyRange*percent, 100*bestSpans/bestSpans(normal), 'Marker', 's');
plot(varyRange*percent, 100*bestChords/bestChords(normal), 'Marker', 'd');
plot(varyRange*percent, 100*bestNP2/bestNP2(normal), 'Marker', '^');
plot(varyRange*percent, 100*bestNP3/bestNP3(normal), 'Marker', 'v');
plot(varyRange*percent, 100*bestTS2/bestTS2(normal), 'Marker', '>');
plot(varyRange*percent, 100*bestTS3/bestTS3(normal), 'Marker', '<');

ylabel('Change in parameter optimal value [%]');
ylim([0 200]);


hold off;

%}