clc;
clear;

%% PLANETOOLS 2022 ~~ RUBBER ~~
planeToolsDirectory = 'D:\Extracurricular\ADT\PlaneTools22\Rubber';

%% Aircraft
airplane = 'Summer21';
plane = ReadPlane(airplane, planeToolsDirectory);

%% Environment parameters
environment.dens = 1.225;

%% Scoring Assumptions
nEggsDroppedBest = 16;
nLapsBest = 16;
nEggsBrokenBest = 1;

bestScore = 0;
planeNo = 0;

%% TRADE STUDY PARAMETERS
chords = .0254*[4 6 8 10 12 14 16 18 20 24];
noPayloads = [1 4 6 8 10 12 14 16 18 20 26 32 48 54];
pitchAngles = (pi/180)*[30 34 38 45];
%pitchAngles = [.5];
Kvs = (2*pi/60)*[180 200 260 300 350 400 500 800];
%Kvs = [52.36];
voltages = 3.7*[2 3 4 6 8 10];
%voltages = [14.8];
diameters = .0254*[6 8 10 12 14 16 18 20 24];

noOfPlanes = length(chords)*length(noPayloads)*length(pitchAngles)*length(Kvs)*length(voltages)*length(diameters);

%% Chord
for chord = chords
    plane.c = chord;

    %% Payloads
    for noPayload = noPayloads
        plane.nPayloads3 = noPayload;

        %% Pitch & Diameter
        for pitchAngle = pitchAngles
            plane.P = pitchAngle;

            for diameter = diameters
                plane.D3 = diameter;

                %% Motor
                for Kv = Kvs
                    plane.Kv = Kv;
                    plane.Rt3 = 2 * ((.000504*plane.Kv) + .019);

                    %% Battery
                    for voltage = voltages
                        
                        if planeNo == 1226
                            a = 1;
                        end
                        plane.Vb = voltage;

                        %% Fuselage, Landing Gear, and Tail Adjustments
                        plane.lFuse = .2 + (.045*plane.nPayloads3);
                        plane.bTail = .4*plane.b;
                        plane.hTail = .2*plane.b;
                        plane.cTail = .1*plane.b;
                        plane.hWing = .2*plane.b;
                        
                        try
                            %% Mass and drag buildups
                            plane.m3 = massBuildup(plane, 3).m;
                            plane.CD03 = dragBuildup(plane, 3).CD0;
                            plane.CD02 = plane.CD03;
                            plane.e = dragBuildup(plane, 3).e;
                            
                            plane.nPayloads2 = 0;
                            plane.m2 = massBuildup(plane, 2).m;
                            
                            %% MISSION 2
                            missionNo = 2;
                            TOthrottle = 1;
                            crzThrottle = 1;
                            turnThrottle = .9;
                            turnCL = .7;

                            V = [];

                            %% THE MISSION
                            % Takeoff
                            segment = FlightModel(plane,environment,missionNo,TOthrottle,0,'takeoff',0);
                            VTO = segment.V;
                            TOFL = segment.TOFL;

                            % Straight
                            segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',VTO);
                            Vcruise = segment.V;
                            Pel = segment.Pel;

                            % Turn
                            segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn180',Vcruise);
                            t180 = segment.t;

                            %% MISSION 3
                            missionNo = 3;
                            TOthrottle = 1;
                            crzThrottle = 1;
                            turnThrottle = .9;
                            turnCL = .7;

                            V = [];

                            %% THE MISSION
                            % Takeoff
                            segment = FlightModel(plane,environment,missionNo,TOthrottle,0,'takeoff',0);
                            VTO = segment.V;
                            TOFL = segment.TOFL;

                            % Straight
                            segment = FlightModel(plane,environment,missionNo,crzThrottle,0,'halfStraightaway',VTO);
                            Vcruise = segment.V;
                            Pel = segment.Pel;

                            % Turn
                            segment = FlightModel(plane,environment,missionNo,turnThrottle,turnCL,'turn180',Vcruise);
                            t180 = segment.t;

                            %% Score mission
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
                            score = (nEggsDropped*nLaps*(1+nEggsBrokenBest))/(nEggsDroppedBest*nLapsBest*(1+nEggsBrokenBest));

                            % Compare to previous best scoring aircraft
                            % If better, save score, no of laps, no of eggs, time, chord, no of
                            % payloads, pitch angle, Kv, voltage, and diameter
                            if score > bestScore
                                bestPlane.no = planeNo;
                                bestPlane.score = score;
                                bestPlane.missionTime = missionTime;
                                bestPlane.nLaps = nLaps;
                                bestPlane.chord = plane.c;
                                bestPlane.nEggsDropped = nEggsDropped;
                                bestPlane.pitchAngle = plane.P;
                                bestPlane.diameter = plane.D3;
                                bestPlane.Kv = plane.Kv;
                                bestPlane.voltage = plane.Vb;
                                bestScore = score;
                            end
                        end
                        planeNo = planeNo + 1
                    end
                end
            end
        end
    end
end
