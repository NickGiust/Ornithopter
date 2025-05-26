function [scores, parameters] = 









%% MASS BUILDUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taken from massBuildup function in PlaneTools22 before the 2022-23 year




%% PARAMETERS

% Fuselage
wFuse = dFuse;
hFuse = dFuse;

% addpath to access materials database
addpath('/Users/jackahrens/Desktop/ADT/PlaneTools22/ComponentLibrary')

fuseMat = 'carbon';
load('material.mat', fuseMat); % figure out how to do this without using readPlane
fuseMaterial = eval(fuseMat);
adFuse = fuseMaterial.density;

pEllipsoid = 1.6075;
aEllipsoid = lFuse/2; bEllipsoid = dFuse/2; cEllipsoid = bEllipsoid;
SwetFuse = 4*pi* ( (aEllipsoid.^pEllipsoid.*bEllipsoid.^pEllipsoid + ...
    aEllipsoid.^pEllipsoid.*cEllipsoid.^pEllipsoid + ...
    bEllipsoid.^pEllipsoid.*cEllipsoid.^pEllipsoid) / 3 ).^(1/pEllipsoid);
% surface area of ellipsoidal fuselage
% from https://en.wikipedia.org/wiki/Ellipsoid as well

% Afuse = 2*(lFuse.*wFuse)+1.5.*(lFuse.*hFuse); % replace with Sfuse for
% ellipsoid approximation
mFuse = adFuse*SwetFuse;

% Wings
wingMat = 'carbon';
load('material.mat', wingMat);
wingMaterial = eval(wingMat);
adWing = wingMaterial.density; % area density?

Awing = 2.2*b.*c;
mWing = adWing*Awing;

Vh = 0.75;
Vv = 0.06;
% bh = 0.4*b;
ARh = 3.5;
lh = 0.6*lFuse;
% Vh = Sh*lh / (S*c)
Sh = Vh.*S.*c./lh;
bh = sqrt(ARh*Sh);
ch = Sh./bh;
lv = lh;
% Vv = Sv*lv / (S*b)
Sv = Vv.*S.*b./lv;
cv = ch;
bv = Sv./cv;
ARv = bv.^2./Sv;

% Tail
cTail = ch;
bTail = bh;
hTail = bv;

tailMat = 'foamfiberglass';
load('material.mat', tailMat);
tailMaterial = eval(tailMat);
adTail = tailMaterial.density;

Atail = 2.2.*cTail.*(bTail+hTail);
mTail = adTail*Atail;

% Motor
motorType = 'Hacker_A50_10L_Turnado_V3';
load('motor.mat', motorType);
motor = eval(motorType);
mMotor = motor.m;

% Battery
% if mission == 2
%     load('battery.mat', plane.batteryType2);
%     battery = eval(plane.batteryType2);
%     mBat = battery.mPerCell*plane.nSeries2*plane.nParallel2;
% elseif mission == 3
%     load('battery.mat', plane.batteryType3);
%     battery = eval(plane.batteryType3);
%     mBat = battery.mPerCell*plane.nSeries3*plane.nParallel3;
% end

batteryType = 'Thunder_LiPo_4400';
nSeries = 6;
nParallel = 1;
load('battery.mat', batteryType);
battery = eval(batteryType);
mBat = battery.mPerCell*nSeries*nParallel;


% Landing gear
lgType = 'bow';
hWing = 0.15;

if strcmp(lgType, 'strut')
    Alg = (hWing^3)*(.15^2);
end

if strcmp(lgType, 'bow')
    Alg = hWing*.2*3.5;
end

lgMat = 'carbon';
load('material.mat', lgMat);
lgMaterial = eval(lgMat);
daLg = lgMaterial.density;

mLg = daLg*Alg;
mLg = (AR.*0 + 1).*mLg;


%% Mass Total
mStructuralComponents = [mFuse mWing mTail mLg];
mStructural = (mFuse + mWing + mTail + mLg) .* buildTechnologyFactor;

mPayload = VWater; % water has 1 kg mass for every liter
mEmpty = mStructural + mBat + mMotor;
mFull = mEmpty + mPayload;











%% TAKEOFF

WFull = mFull*g; % in N
WEmpty = mEmpty*g; % in N

VstallFull = sqrt(2.*WFull ./ (rho.*S.*CLmax)); % set liftoff velocity, also to check cruise condition 
VstallEmpty = sqrt(2.*WEmpty ./ (rho.*S.*CLmax)); % to check cruise condition


VLO = 1.2*VstallFull; % liftoff velocity after takeoff roll is a margin above stall

phi = (16*hWing./b).^2 ./ (1 + (16*hWing./b).^2);

CLtakeoffRun = 0.7.*CLmax; % average lift coefficient during takeoff run

k = 1 / (pi.*e.*AR);

qTakeoffRun = 0.5*dens.*(0.7*VLO).^2; % average velocity during takeoff run

dragTakeoff = (CD0 + phi.*k.*CLtakeoffRun.^2) .* qTakeoffRun .* S;

liftTakeoff = CLtakeoffRun .* qTakeoffRun .* S; % average lift force during takeoff


% compute TOFL
TOFL = 1.44.*WFull.^2 ./ (dens*g.*S.*CLmax .* (TStatic - (dragTakeoff + Crr.*(WFull - liftTakeoff))));
% see 261 final exam equation sheet


% compute TOFL energy
energyTakeoff = dragTakeoff.*TOFL; % energy consumed is distance*drag

clear VLO phi qTakeoffRun dragTakeoff liftTakeoff


%% Laps
% NOTE: ENERGY IS ADDED UP AS FORCE*DISTANCE, THEN EFFICIENCY FACTORS ARE
% APPLIED AT THE END

%{
bList = bLimit; % span will always be the maximum to reduce induced drag
ARList = 4:2:8; % aspect ratio between 3 and 15, seems reasonable
CLmaxList = 2:0.1:2.2; % maximum lift coefficient between 1.4 and 2.2
CLCruiseFullList = 0.4:0.2:0.8; % average cruise lift coefficient while heavy
CLCruiseEmptyList = 0.2:0.2:0.6; % average cruise lift coefficient after dumping water
TStaticList = 40:40:120; % static thrust, N ( = 1 kg to 20 kg)
TDynamicShareFullList = 0.1:0.1:0.5; % drag in cruise while heavy
TDynamicShareEmptyList = 0.1:0.1:0.5; % drag in cruise after dumping water
nTurnFullList = 2:2:6; % number of G's pulled in turns while full, assuming level
nTurnEmptyList = 4:2:10; % number of G's pulled in turns while empty, assuming level
VWaterList = 1:2:5; % volume of water payload, Liters
numLapsList = 7:1:12; % total number of laps flown
numLapsFull = ceil(numLapsList/2); % at least half of laps must be with water
numLapsEmpty = numLapsList - numLapsFull; % at max half of laps are empty
buildTechnologyFactorList = 0.9:0.1:1.1; % multiplier on aircraft empty weight

% Is having both CLCruise and TDynamicShare redundant?
% Since T = D = 0.5*rho*U^2*S*(CL^2/(pi*e*AR) + CDO)
% density is fixed, speed is unknown, S is swept over, CL is unknown
% e estimated as a constant value, AR is swept over, and CDO = f(geometry)
% which is known (assuming speed is estimated for Re calculation)
% In summary, thrust/drag, speed, and lift coefficient are unknown.
% Therefore, two of thrust/speed/lift must be set to have a determinate
% (steady state) cruise condition
% Choose to set lift coefficient

%}


%% THRUST/DRAG

% TDynamicFull = TDynamicShareFull.*TStatic;
% DynamicEmpty = TDynamicShareEmpty.*TStatic;


%% STRAIGHTS

speedStraightFull = sqrt(WFull ./ (0.5 * rho .* S .* CLCruiseFull));
speedStraightEmpty = sqrt(WEmpty ./ (0.5 * rho .* S .* CLCruiseEmpty));

CDStraightFull = CLCruiseFull.^2.*k + CD0;
CDStraightEmpty = CLCruiseEmpty.^2.*k + CD0;

DStraightFull = CDStraightFull .* S * 0.5 * rho .* speedStraightFull.^2;
DStraightEmpty = CDStraightEmpty .* S * 0.5 * rho .* speedStraightEmpty.^2;


%% TURNS

turnSpeedFraction = 0.85;
% take average speed in turns to be 85% of speed in cruise due to
% deceleration

speedTurnFull = turnSpeedFraction*speedStraightFull;
speedTurnEmpty = turnSpeedFraction*speedStraightEmpty;

qTurnFull = 0.5*rho*speedTurnFull.^2;
qTurnEmpty = 0.5*rho*speedTurnEmpty.^2;

CLTurnFull1G = WFull./(S.*qTurnFull);
CLTurnEmpty1G = WEmpty./(S.*qTurnEmpty);

CLTurnFull = CLTurnFull1G.*nTurnFull;
CLTurnEmpty = CLTurnEmpty1G.*nTurnEmpty;

CDTurnFull = CLTurnFull.^2.*k + CD0;
CDTurnEmpty = CLTurnEmpty.^2.*k + CD0;

DTurnFull = CDTurnFull .* S .* qTurnFull;
DTurnEmpty = CDTurnEmpty .* S .* qTurnEmpty;

rTurnFull = 2*qTurnFull ./ (rho*g*sqrt(nTurnFull.^2-1));
rTurnEmpty = 2*qTurnEmpty ./ (rho*g*sqrt(nTurnEmpty.^2-1));

clear qTurnFull qTurnEmpty CLTurnFull1G CLTurnEmpty1G

%% DISTANCE and TIME

% circumference = pi*d, half circumference = pi*d/2 = pi*r
distance180Full = pi*rTurnFull;
distance180Empty = pi*rTurnEmpty;

% distance = rate * time --> time = distance/rate

timeHalfStraightFull = distanceHalfStraight ./ speedStraightFull;
timeHalfStraightEmpty = distanceHalfStraight ./ speedStraightEmpty;

time180Full = distance180Full ./ speedTurnFull;
time180Empty = distance180Empty ./ speedTurnEmpty;


%% ENERGY and TIME CALCULATION

energyHalfStraightFull = DStraightFull .* distanceHalfStraight;
energyHalfStraightEmpty = DStraightEmpty .* distanceHalfStraight;

energy180Full = DTurnFull .* distance180Full;
energy180Empty = DTurnEmpty .* distance180Empty;

energyTurnReactionTimeFull = speedStraightFull .* turnReactionTime .* DStraightFull;
energyTurnReactionTimeEmpty = speedStraightEmpty .* turnReactionTime .* DStraightFull;

energyLapFull = 4*energyHalfStraightFull + 4*energy180Full + 2*energyTurnReactionTimeFull;
energyLapEmpty = 4*energyHalfStraightEmpty + 4*energy180Empty + 2*energyTurnReactionTimeEmpty;

% percent of each lap's energy dedicated to each segment of the lap
energyShareFullStraights = 4*energyHalfStraightFull ./ energyLapFull;
energyShareFullTurns = 4*energy180Full ./ energyLapFull;
energyShareFullTurnReactionTime = 2*energyTurnReactionTimeFull ./ energyLapFull;

energyShareEmptyStraights = 4*energyHalfStraightEmpty ./ energyLapEmpty;
energyShareEmptyTurns = 4*energy180Empty ./ energyLapEmpty;
energyShareEmptyTurnReactionTime = 2*energyTurnReactionTimeEmpty ./ energyLapEmpty;


timeLapFull = 4*timeHalfStraightFull + 4*time180Full;
timeLapEmpty = 4*timeHalfStraightEmpty + 4*time180Empty;


energyFull = numLapsFull .* energyLapFull;
energyEmpty = numLapsEmpty .* energyLapEmpty;

energyTotalMechanical = energyTakeoff + energyFull + energyEmpty;

energyElectric = energyTotalMechanical / etaTotal;


timeFull = numLapsFull .* timeLapFull;
timeEmpty = numLapsEmpty .* timeLapEmpty;

timeTotal = timeFull + timeEmpty;





%% TRIM RESULTS TO TOFL, ENERGY/EFFICIENCY, TIME, AND STALL

score = VWater .* numLaps; % water weight = water volume between kg and L

% set scores to 0 if TOFL not met
score(TOFL > TOFLLimit) = 0;

% set scores to 0 if battery energy and SoC exceeded
score(energyElectric > (1 - batteryEndSoC) * energyLimitJ) = 0;

% set scores to 0 if the total lap time exceeds the 10 minute time limit
score(timeTotal > timeLimit) = 0;

% set scores to 0 if the straightaway cruise velocity is not at least
% somewhat higher than the stall velocity
score(speedStraightFull < 1.2 * VstallFull) = 0;
score(speedStraightEmpty < 1.2 * VstallEmpty) = 0;

% set scores to 0 if the turns require a higher lift coefficient than CLmax
score(CLTurnFull > 0.9*CLmax) = 0;
score(CLTurnEmpty > 0.9*CLmax) = 0;


%% OUTPUT RESULTS

clc;

[bestScore, bestScoreIndex] = max(score, [], 'all', 'linear');
% bestScoreIndex reported as a single value
fprintf('\nThe best score is %.2f out of %.1f million total aircraft combinations\n\n', bestScore, numPlanes/(1e6))
fprintf('It completes the mission in %.1f seconds, or %.2f minutes, by flying %g laps\n', timeTotal(bestScoreIndex), timeTotal(bestScoreIndex)/60, numLaps(bestScoreIndex))
fprintf('The water volume for best score is %.2f L. Range was %.2f to %.2f L\n', VWater(bestScoreIndex), VWaterList(1), VWaterList(end))
fprintf('The number of laps for best score is %.f. Range was %.f to %.f laps\n\n', numLaps(bestScoreIndex), numLapsList(1), numLapsList(end))


fprintf('This aircraft has an empty (non-payload weight) of %.2f kg\n', mEmpty(bestScoreIndex))
fprintf('The component weights are:\n')
fprintf('           Wing: %1.2f kg\n', mTail(bestScoreIndex))
fprintf('          Tails: %1.2f kg\n', mWing(bestScoreIndex))
fprintf('       Fuselage: %1.2f kg\n', mFuse(bestScoreIndex))
fprintf('   Landing Gear: %1.2f kg\n', mLg(bestScoreIndex))
fprintf('        Battery: %1.2f kg\n', mBat)
fprintf('          Motor: %1.2f kg\n', mMotor)

fprintf('\nThe wingspan was set to the maximum of 4 ft = %.2f m to minimize induced drag\n', 4*0.3048)
fprintf('The AR for best score is %.2f. Range was %.2f to %.2f\n', AR(bestScoreIndex), ARList(1), ARList(end))
fprintf('Therefore, the wing area is %.3f m^2 and the chord is %.2f m\n', S(bestScoreIndex), c(bestScoreIndex))
fprintf('Additionally, the fuselage is VERY roughly %.2f m long with a diameter of %.2f m\n', lFuse(bestScoreIndex), dFuse(bestScoreIndex))
fprintf('This geometry has a parasite drag coefficient of %.3f\n\n', CD0(bestScoreIndex))

fprintf('The CLmax for best score is %.1f. Range was %.1f to %.1f\n', CLmax(bestScoreIndex), CLmaxList(1), CLmaxList(end))
fprintf('The static thrust for best score is %g N or %.2f kg. Range was %.f to %.f N\n', TStatic(bestScoreIndex), TStatic(bestScoreIndex)/g, TStaticList(1), TStaticList(end))
fprintf('With these parameters, the takeoff field length is %.2f m (limit = %.2f m), and stall velocity is %.1f m/s\n\n', TOFL(bestScoreIndex), TOFLLimit, VstallFull(bestScoreIndex))


fprintf('The CLCruiseFull for best score is %.2f. Range was %.2f to %.2f\n', CLCruiseFull(bestScoreIndex), CLCruiseFullList(1), CLCruiseFullList(end))
fprintf('This gives speeds w/ water of %.1f m/s in cruise and %.1f m/s in turns\n', speedStraightFull(bestScoreIndex), speedTurnFull(bestScoreIndex))
fprintf('The CLCruiseEmpty for best score is %.2f. Range was %.2f to %.2f\n', CLCruiseEmpty(bestScoreIndex), CLCruiseEmptyList(1), CLCruiseEmptyList(end))
fprintf('This gives speeds w/o water of %.1f m/s in cruise and %.1f m/s in turns\n\n', speedStraightEmpty(bestScoreIndex), speedTurnEmpty(bestScoreIndex))


fprintf('The turn G''s w/ water for best score was %.1f. Range was %.1f to %.1f\n', nTurnFull(bestScoreIndex), nTurnFullList(1), nTurnFullList(end))
fprintf('The turn radius w/ water for best score was %.1f m, pulling a CL of %.2f (CLmax = %.2f)\n', rTurnFull(bestScoreIndex), CLTurnFull(bestScoreIndex), CLmax(bestScoreIndex))
fprintf('The turn G''s w/o water for best score was %.1f. Range was %.1f to %.1f\n', nTurnEmpty(bestScoreIndex), nTurnEmptyList(1), nTurnEmptyList(end))
fprintf('The turn radius w/o water for best score was %.1f m, pulling a CL of %.2f (CLmax = %.2f)\n', rTurnEmpty(bestScoreIndex), CLTurnEmpty(bestScoreIndex), CLmax(bestScoreIndex))


%% DEBUGGING
% figure out why tails weigh nearly as much as (or more than!) the wing



