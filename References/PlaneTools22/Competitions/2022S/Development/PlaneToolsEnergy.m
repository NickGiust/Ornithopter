clc;
clear;
warning('off','all');

%% PLANETOOLS 2022 ~~ RUBBER 2.0 ~~
planeToolsDirectory = 'C:\Users\jacka\OneDrive\Documents\GitHub\PlaneTools22';

%% Aircraft
seedPlane = 'ScoreAnalysis22';
pl = ReadPlane(seedPlane, '22', planeToolsDirectory);

%% Environment parameters
dens = 1.225;

%% Estimated parameters ("fudge" factors)
% Proportion of static thrust produced at max efficiency
dynamicShare = .35;

% Value during turns
dynamicShareTurn = .45;

% Total propulsive efficiency
eta = .65;

% Maximum CL on takeoff
CLmax = 2.2;

% Maximum CL in turns
turnCL = 2.1;

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
tLoadSyringe = .2;
tUnloadSyringe = .4;
tLoadPackage = 3;

% Technological advancement (total ability to reduce weight and drag)
tech = 1;

TOFLlimit = 7.62;

%% Scoring Assumptions
syringeSpeedBest = 10;
nLapsBest = 15;
tGMbest = 10;

%% TRADE STUDY PARAMETERS
chords = 0.25:0.05:0.6;
%chords = .25;
noPackages = 8:1:12;
%noPackages = 10;
noSyringes = 400:25:600;
%noSyringes = 100;
staticThrust3 = 40:10:90;
%staticThrust3 = 60;
staticThrust2 = 180:20:260;
%staticThrust2 = 240;
%spans = 1.8:0.2:2.4;
spans = 2.4;
wingletHeightsRelative = 0:0.05:0.35;
%wingletHeightsRelative = 0;

totalNumberPlanes = length(chords)*length(noPackages)*length(noSyringes)*...
    length(staticThrust3)*length(staticThrust2)*length(spans)*length(wingletHeightsRelative);

fprintf('%d total planes\n\n', totalNumberPlanes)
fprintf('start -------------------- end\n')
fprintf('      ')
asteriskCounter = 0;
asteriskCounterPrevious = 0;

scores = [];
noFailure = 1;
bestScore = 0;
planeNo = 0;
noScoringPlanes = 0;

varyRange = 10;
variable = 'Best syringe speed [syringes/s]';

%% Score sensitivity parameter
for scoreParam = varyRange
    syringeSpeedBest = scoreParam;
    
    bestScore = 0;
    %planeNo = 0;
    
    %% Chord
    for wingspan = spans
        pl.b = wingspan;

        for wingletHeightRelative = wingletHeightsRelative
            
            wingletHeight = wingletHeightRelative*wingspan;
            pl.hwinglet = wingletHeight;
        
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
       
                            planeNo = planeNo + 1;
                            planesFraction = planeNo/totalNumberPlanes;
                            asteriskCounterPrevious = asteriskCounter;
                            asteriskCounter = floor(planesFraction*100/5);
                            if asteriskCounter > asteriskCounterPrevious
                                fprintf('*')
                            end
                            
                            %% Fuselage, Landing Gear, and Tail Adjustments
                            S = pl.b*pl.c;
                            Cvt = .08;
                            Cht = .08;
                            pl.lFuse = max([(pl.b/2) (.3+(.16*pl.nPayloads3))]);
                            Svt = Cvt*pl.b/(pl.lFuse*.7);
                            Sht = Cht*pl.b/(pl.lFuse*.7);
                            pl.cTail = (Sht/2.5)^.5;
                            pl.bTail = 2.5*pl.cTail;
                            pl.hTail = Svt/pl.cTail;
                            pl.hWing = .2*pl.b;

                            if pl.lFuse > 2.44
                                pl.lFuse = 2.44;
                                pl.hFuse = .14 * (pl.nPayloads3/12);
                            end

                            %% Mass and drag buildups
                            m = massBuildup(pl, 2).m*tech;
                            CD0 = dragBuildup(pl, 2).CD0*tech;
                            e = dragBuildup(pl, 3).e;
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
                            Vturn = ((2*dynamicShareTurn*Tstatic2)/(dens*S*(CD0+(k*(turnCL^2)))))^.5;
                            Lturn = 0.5*dens*(Vturn^2)*S*turnCL;
                            loadFactor = Lturn/W;
                            S = S/1.2;
                            
                            if loadFactor > pl.nStruct
                                loadFactor = pl.nStruct;
                            end
                            
                            r = (Vturn^2)/(9.81*loadFactor);
                            d180 = pi*r;
                            t180 = d180/Vturn;
                            
                            loadFactor2 = loadFactor;

                            %% Takeoff (distance)
                            m2 = m;
                            Vstall = sqrt((2*W)/(CLmax*dens*S));
                            VTO = 1.2*Vstall;
                            gfx = ((16*pl.hWing/pl.b)^2)/(1+((16*pl.hWing/pl.b)^2));

                            % Acceleration = force/mass
                            staticThr = Tstatic2/m;
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
                            % Takeoff field length
                            if isempty(x) || x(end) > TOFLlimit
                                continue;
                            else
                                TOFL2 = x(end);
                            end

                            %% Score
                            % Time to fly a lap
                            tLap = (612/Vcruise) + (t180*4);

                            % Add wind drift
                            tWindDrift = (tLap*headwind)/Vcruise;
                            tLap = tLap + tWindDrift;

                            tM2 = 3*tLap;

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
                            m = massBuildup(pl, 3).m*tech;
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
                            Vturn = ((2*dynamicShareTurn*Tstatic3)/(dens*S*(CD0+(k*(turnCL^2)))))^.5;
                            Lturn = 0.5*dens*(Vturn^2)*S*turnCL;
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
                            tGM = (tLoadSyringe + tUnloadSyringe)*pl.nPayloads2 + tLoadPackage*pl.nPayloads3;
                            if tGM > 300
                                continue;
                            end
                            scoreGM = tGMbest/tGM;

                            %% Score plane and compare
                            score = (score1 + score2 + score3 + scoreGM)*noFailure;

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
                                bestPlane.CD0 = CD0;
                                bestPlane.nSyringes = pl.nPayloads2;
                                bestPlane.nPackages = pl.nPayloads3;
                                bestPlane.Tstatic2 = Tstatic2;
                                bestPlane.Tstatic3 = Tstatic3;
                                bestPlane.hWinglet = wingletHeight;
                                bestPlane.hWingletRelative = wingletHeightRelative;
                                
                                bestPlane.syringeSpeed = pl.nPayloads2/tM2;
                                bestPlane.nLaps3 = nLaps;
                                bestPlane.tGM = tGM;

                                bestPlane.tM2 = tM2;
                                bestPlane.TOFL2 = TOFL2;
                                bestPlane.m2 = m2;
                                bestPlane.PF2 = (pl.mPayload2*pl.nPayloads2)/m2;
                                bestPlane.Vcruise2 = Vcruise2;
                                bestPlane.loadFactor2 = loadFactor2;
                                bestPlane.Pel2 = Pel2;

                                bestPlane.tM3 = missionTime;
                                bestPlane.TOFL3 = TOFL3;
                                bestPlane.mTO3 = m;
                                bestPlane.PF3 = (pl.mPayload3*pl.nPayloads3)/m;
                                bestPlane.Vcruise3 = Vcruise3;
                                bestPlane.loadFactor3 = loadFactor;
                                bestPlane.Pel3 = Pel3;

                                bestScore = score;
                            end                   
                            noScoringPlanes = noScoringPlanes + 1;
                            end
                        end
                    end
                end
            end
        end
    scores = [scores bestScore];
    
    end
end

fprintf('\n%d planes scored out of %d total planes\n', noScoringPlanes, totalNumberPlanes)
bestPlane

plot(varyRange, scores);
ylabel('DBF Score');
ylim([3 7]);
xlabel(variable);
