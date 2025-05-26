function [mWiring] = calcMass_Wiring(plane)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INPUTS
b = plane.b;
c = plane.c;

prop2DiamInches = plane.D2; % input
prop3DiamInches = plane.D3;
%prop3DiamInches = str2num(plane.propType3(strfind(plane.propType3,'_')+1:strfind(plane.propType3,'x')-1)); % input

propDiamInches = max([prop2DiamInches prop3DiamInches]);

fuseLength = plane.lFuse;
fuseWidth = plane.wFuse;

%Motor Count
motorCount = plane.nMotors;

%Motor Location, "n", nose, or "w", wing
if motorCount == 1
    motorLoc = "n";
else
    motorLoc = "w";
end

%Continuous Current HERE

%% CALC

fracFusetoWingEntry = 0.25;
fuseToWingEntry = fracFusetoWingEntry * fuseLength;

%Arming Plug
armPlugFuseFrac = 0.6/2.3; %measured VAX
armPlugDistNose = armPlugFuseFrac * fuseLength;

%If wing mount, percent of half span
if motorLoc == "w"
    %motorHalfSpanPerc = 0.3;
    %distance from middle of fuselage to motor on either side
    %take into account: prop radius, fuse width, and propFuseClearancem tip/fuselage
    %clearance
    propFuseClearance = 0.15;
    motorLateralDist = (((propDiamInches/2)/39.37)) + fuseWidth/2 + propFuseClearance;
end

%ESC Count
escCount = motorCount;

%Flap Servos Half Span Percentage
flapServosHalfSpanPerc = 0.2;

%Ailerons Servos Half Span Percentage
ailerServosHalfSpanPerc = 0.6;

servoMassPerLen = 0.0216/1.77;

fracMotorChord = 0.5;

lengthMargin = 1.05;

%--------------------------------------------------------------------

%%%IMPORT GAUGE WEIGHT TABLE
GaugesWeightTable = readtable('Gauges Weight Table.xlsx', 'Sheet', 'Sheet2');

%Weight/Length data: https://en.wikibooks.org/wiki/Engineering_Tables/Standard_Wire_Gauge\
%Outer Diam: https://www.panduit.com/content/dam/panduit/en/products/media/4/54/254/3254/13254.pdf
%Inner Diam: https://www.meridiancableassemblies.com/2021/04/wire-gauge-size-guide/

%-------------------------------------------------------------------
%%%CALCULATING GAUGES
%Gauges
%TO BE IMPLEMENTED, HARD CODED BELOW (IDEALLY FUNCTION OF MAXCURRENT)
motorWireGauge = 8;
servosWireGauge = 26;
batteryWireGauge = motorWireGauge;

%%%CREATING LENGTHS ARRAYS
%A length array holds every individual calculated length
%Empty Lengths Arrays

motorWireLen = [];
servosWireLen = [];
batteryWireLen = [];

%-------------------------------------------------------------------

%%%MOTOR WIRING (ESC-MOTOR) CALCULATIONS
%Motor
if motorLoc == "n" %Nose mount, assume 30% of fuse
    motorWireLen = [motorWireLen 3*(0.3*fuseLength)];
    %multiply by three because three leads
    %assume 30%, talked to nick, need a little to allow CG adj.
elseif motorLoc == "w" %Wing mount
    %Append to length array the motor lateral distance (from fuse to motor)
    %for 3 leads
    motorWireLen = [motorWireLen 3*(motorCount*motorLateralDist)];

    %Append to length array the motor horizontal distance (from somewhere
    %in the chord to the motor position) for 3 leads
    motorWireLen = [motorWireLen  3*(motorCount*fracMotorChord*c)];

    %Append to length array the distance from the electronics to the wing
    %entry point for 3 leads
    motorWireLen = [motorWireLen 3*(fuseToWingEntry)];
else
    error("Invalid motor location input")
end

%%%SERVOS WIRING (ESC-SERVOS) CALCULATIONS
%Servos

%Tail and rudder
servosWireLen = [servosWireLen 3*(2*(0.5*fuseLength))];
%multiply by three because three leads

%Ailerons
servosWireLen = [servosWireLen 3*(2*(ailerServosHalfSpanPerc*b/2))];

%Flaps
servosWireLen = [servosWireLen 3*(2*(flapServosHalfSpanPerc*b/2))];

%%%BATTERY WIRING (BATTERY-ESC) CALCULATIONS
%%BATTERY ESC
%two leads, 7cm each lead
batteryWireLen = [batteryWireLen 2*(escCount*(0.07))];


%%ARMING PLUG
%2 leads (x2), goes and comes back (x2)
batteryWireLen = [batteryWireLen 2*(2*armPlugDistNose)];

%----------------------------------------------------------------------

%%%WEIGHT ESTIMATION: WIRE
motorWireWeight = lengthMargin*sum(table2array(GaugesWeightTable(motorWireGauge, 2)) .* motorWireLen);
servosWireWeight = lengthMargin*sum(servoMassPerLen .* servosWireLen);
batteryWireWeight = lengthMargin*sum(table2array(GaugesWeightTable(batteryWireGauge, 2)) .* batteryWireLen);

%%%WEIGHT ESTIMATION: INSULATION
%Only doing for motor and ESC for now, no reliable values for servos (also, quite small, so somwhat negligible)
%1500 kg/m^3, density for silicone rubber insulator
insulatorDensity = 1500;

%0.001 to convert mm to m
outerRadiusSquared = ((table2array(GaugesWeightTable(motorWireGauge, 4))/2)*0.001)^2;
innerRadiusSquared = ((table2array(GaugesWeightTable(motorWireGauge, 3))/2)*0.001)^2;
motorInsulationCrossSection = pi*(outerRadiusSquared - innerRadiusSquared);

motorInsulationWeight = lengthMargin*(insulatorDensity*(motorInsulationCrossSection * sum(motorWireLen)));

batteryInsulationCrossSection = motorInsulationCrossSection;
batteryInsulationWeight = lengthMargin*(insulatorDensity*(batteryInsulationCrossSection * sum(batteryWireLen)));

%%%WEIGHT ESTIMATION: TOTAL
mWiring = motorWireWeight + servosWireWeight + batteryWireWeight + motorInsulationWeight + batteryInsulationWeight;

end

