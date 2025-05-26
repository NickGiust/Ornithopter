clc;
clear;

%% PLANETOOLS 2022 ~~ RUBBER 2.0 ~~
planeToolsDirectory = 'D:\ADT\PlaneTools22\';

%% Aircraft
airplane = 'Summer21';
plane = ReadPlane(airplane, planeToolsDirectory);

%% Environment parameters
dens = 1.225;

%% Propulsion performance parameters
%% REMEMBER DESIGN SPACE EDGES: TOFL, STRUCTURE DURING TURNS, ASSUMPTIONS
% Proportion of static thrust produced at max efficiency
dynamicShare = .35;
% Value during turns
dynamicShareTurn = .5;
% Takeoff EFFECTIVE (unsteady) dynamic share
dynamicShareTO = .8;
% Total propulsive efficiency
eta = .6;
% Structural load factor limit
nStruct = 5;
% Maximum CL on takeoff
CLmax = 2.1;
% Maximum CL in turns
turnCL = 1.5;
% Headwind
headwind = 1;

%% Basic parameters
b = plane.b;
hWing = plane.hWing;
CLgroundRoll = plane.CLgroundRoll;
rr = plane.rr;

%% Scoring Assumptions
TOFLlimit = 9.14;
nEggsDroppedBest = 20;
nLapsBest = 20;
nEggsBrokenBest = 1;
tM2Best = 50;

bestScore = 0;
planeNo = 0;

%% TRADE STUDY PARAMETERS
chords = .0254*[4 5 6 7 8 9 10 11 12];
%chords = .0254*12;
noPayloadsHalf = [2 3 4 5 6 7 8 9 10 11];
%noPayloads = 3;
staticThrust = [5 10 20 30 40 50 65 75];
%staticThrust = 67;

%% Chord
for chord = chords
    plane.c = chord;

    %% Payloads
    for noPayloadHalf = noPayloadsHalf
        plane.nPayloads3 = noPayloadHalf;
        
        for Tstatic2 = staticThrust
            
            for Tstatic3 = staticThrust
     
                %% Fuselage, Landing Gear, and Tail Adjustments
                plane.lFuse = .2 + (.045*plane.nPayloads3);
                plane.bTail = .4*plane.b;
                plane.hTail = .3*plane.b;
                plane.cTail = .1*plane.b;
                plane.hWing = .2*plane.b;

                try
                    %% Mass and drag buildups
                    W = (massBuildup(plane, 2).m)*9.81;
                    CD0 = dragBuildup(plane, 2).CD0;
                    e = dragBuildup(plane, 3).e;

                    S = plane.b*plane.c;
                    AR = (plane.b^2)/S;
                    k = 1/(pi*e*AR);
                    S = (plane.b-plane.wFuse)*plane.c;

                    %% MISSION 2
                    % Straight
                    Vcruise = (((dynamicShare*Tstatic2)+((((dynamicShare^2)*(Tstatic2^2))-(4*k*(W^2)*CD0))^.5))/(CD0*dens*S))^.5;
                    Pel = (dynamicShare*Tstatic2*Vcruise)/eta;
                    Vcruise2 = Vcruise;
                    
                    if ~isreal(Vcruise)
                        error('Nonreal cruise condition.');
                    end

                    % Turn
                    Vturn = ((2*dynamicShareTurn*Tstatic2)/(dens*S*(CD0+(k*(turnCL^2)))))^.5;
                    Lturn = 0.5*dens*(Vturn^2)*S*turnCL;
                    loadFactor = Lturn/W;
                    r = (Vturn^2)/(9.81*loadFactor);
                    d180 = pi*r;
                    t180 = d180/Vturn;

                    if loadFactor > nStruct
                        error('Structure warning: aircraft structural load factor limit exceeded during turn.');
                    end

                    % Takeoff
                    mTO = W/9.81;
                    Vstall = sqrt((2*W)/(CLmax*dens*S));
                    VTO = 1.2*Vstall;
                    gfx = ((16*hWing/b)^2)/(1+((16*hWing/b)^2));

                    % Time to takeoff
                    tspan = 0:.1:30;
                    A = (dynamicShareTO*Tstatic2)/mTO;
                    B = (dens*S*(CD0 + (k*gfx*(CLgroundRoll^2)) - (rr*CLgroundRoll)))/(2*mTO);
                    C = (rr*W)/mTO;
                    [t,V] = ode45((@(t,V) A-C-(B*(V^2))), tspan, headwind);

                    TOindex = find( V > VTO, 1 );
                    tTO = t(TOindex);

                    % V as a function of t up to takeoff speed
                    tspan = 0:.1:tTO;
                    [t,V] = ode45((@(t,V) A-C-(B*(V^2))), tspan, headwind);

                    % Distance as a function of t up to takeoff speed
                    d = cumsum(diff(t).*(V(1:end-1)-headwind));

                    % Takeoff field length in FEET
                    TOFL = d(end);

                    if TOFL > TOFLlimit
                        error('Takeoff field length too long.');
                    end

                    %% Score
                    % Time to fly a lap
                    tLap = (612/Vcruise) + (t180*4);

                    tM2 = 3*tLap;
                    
                    if planeNo == 5420
                        a = 1;
                    end
                    if Pel > 2300
                        error('Electric power limit exceeded');
                    end

                    % Energy to fly a lap
                    Elap = Pel*tLap;

                    % No of laps to 100 Wh
                    nLapsEnergy = floor(100*3600/Elap);

                    % Score 0 if less than 3 laps are completed
                    if nLapsEnergy < 3
                        score2 = 0;
                    else
                        score2 = 1*(tM2Best/tM2);
                    end

                    %% MISSION 3
                    W = (massBuildup(plane, 3).m)*9.81;
                    CD0 = dragBuildup(plane, 3).CD0;

                    % Straight
                    Vcruise = (((dynamicShare*Tstatic3)+((((dynamicShare^2)*(Tstatic3^2))-(4*k*(W^2)*CD0))^.5))/(CD0*dens*S))^.5;
                    Pel = (dynamicShare*Tstatic3*Vcruise)/eta;
                    Vcruise3 = Vcruise;
                    
                    if ~isreal(Vcruise)
                        error('Nonreal cruise condition.');
                    end
                    
                    % Turn
                    Vturn = ((2*dynamicShareTurn*Tstatic3)/(dens*S*(CD0+(k*(turnCL^2)))))^.5;
                    Lturn = 0.5*dens*(Vturn^2)*S*turnCL;
                    loadFactor = Lturn/W;
                    r = (Vturn^2)/(9.81*loadFactor);
                    d180 = pi*r;
                    t180 = d180/Vturn;

                    if loadFactor > nStruct
                        error('Structure warning: aircraft structural load factor limit exceeded during turn.');
                    end

                    % Takeoff
                    mTO = W/9.81 + (plane.mPayload3*plane.nPayloads3);
                    Vstall = sqrt((2*W)/(CLmax*dens*S));
                    VTO = 1.2*Vstall;
                    gfx = ((16*hWing/b)^2)/(1+((16*hWing/b)^2));

                    % Time to takeoff
                    tspan = 0:.1:30;
                    A = (dynamicShareTO*Tstatic3)/mTO;
                    B = (dens*S*(CD0 + (k*gfx*(CLgroundRoll^2)) - (rr*CLgroundRoll)))/(2*mTO);
                    C = (rr*W)/mTO;
                    [t,V] = ode45((@(t,V) A-C-(B*(V^2))), tspan, headwind);

                    TOindex = find( V > VTO, 1 );
                    tTO = t(TOindex);

                    % V as a function of t up to takeoff speed
                    tspan = 0:.1:tTO;
                    [t,V] = ode45((@(t,V) A-C-(B*(V^2))), tspan, headwind);

                    % Distance as a function of t up to takeoff speed
                    d = cumsum(diff(t).*(V(1:end-1)-headwind));

                    % Takeoff field length in FEET
                    TOFL3 = d(end);

                    if TOFL3 > TOFLlimit
                        error('Takeoff field length too long.');
                    end

                    %% Score
                    % Time to fly a lap
                    tLap = (612/Vcruise) + (t180*4);

                    % No of laps to 10 min
                    nLapsTime = floor(600/tLap);

                    % Energy to fly a lap
                    Elap = Pel*tLap;

                    % No of laps to 100 Wh
                    nLapsEnergy = floor(100*3600/Elap);

                    % Whichever is less is the number of laps
                    nLaps = min([nLapsTime nLapsEnergy]);

                    if nLaps < 0
                        error('Invalid number of laps.');
                    elseif nLaps < nLapsTime
                        missionTime = tLap*nLaps;
                    else
                        missionTime = 600;
                    end

                    % If no of eggs carried is greater than no of laps, eggs dropped is number
                    % of laps
                    if 2*plane.nPayloads3 < nLaps
                        nEggsDropped = 2*plane.nPayloads3;
                    else
                        nEggsDropped = nLaps;
                    end

                    % No of eggs broken
                    nEggsBroken = ceil(nEggsDropped/5);

                    % Score mission
                    score3 = (nEggsDropped*nLaps*(1+nEggsBrokenBest))/(nEggsDroppedBest*nLapsBest*(1+nEggsBrokenBest));

                    %% Score plane and compare
                    score = score2 + score3;

                    % Compare to previous best scoring aircraft
                    % If better, save score, no of laps, no of eggs, time, chord, no of
                    % payloads, pitch angle, Kv, voltage, and diameter
                    if score > bestScore
                        bestPlane.no = planeNo;
                        bestPlane.score2 = score2;
                        bestPlane.score3 = score3;
                        bestPlane.score = score;
                        bestPlane.tM2 = tM2;
                        bestPlane.TOFL2 = TOFL;
                        bestPlane.tM3 = missionTime;
                        bestPlane.nLaps3 = nLaps;
                        bestPlane.TOFL3 = TOFL3;
                        bestPlane.mTO3 = mTO;
                        bestPlane.Tstatic2 = Tstatic2;
                        bestPlane.Vcruise2 = Vcruise2;
                        bestPlane.Tstatic3 = Tstatic3;
                        bestPlane.Vcruise3 = Vcruise3;
                        bestPlane.chord = plane.c;
                        bestPlane.nEggsDropped = nEggsDropped;
                        bestPlane.nEggsBroken = nEggsBroken;
                        bestScore = score;
                    end
                end
                planeNo = planeNo + 1
            end
        end
    end
end

