function [tubePlyArea, shearWebVolume, plane] = sparSizerGM_2(plane)

% Jack Ahrens - jahrens@usc.edu - 10/16/22

plane.validTubes = true;

% %{

b = plane.b;
% AR = plane.AR;
tcRatio = plane.tcRatio;
GMloadKG = plane.GMloadingMultiplier * plane.seed_m2; % mass, in kg, put on plane during GM

% inputs
semispan = b / 2; % [in] half-span (single panel)

ds = 0.01; % numerical increment, m
x = (0:ds:semispan);
numPoints = length(x);
% buttLine = x/semispan; % fraction of way between root and tip

% cavg = 0.5*(croot + ctip)
% taperRatio = ctip/croot
% ctip = croot*taperRatio
% cavg = 0.5*(croot + croot*taperRatio)
% cavg = 0.5*croot*(1+taperRatio)
% croot = 2*cavg / (1 + taperRatio)
% croot = 2*cavg / (1 + taperRatio);
% ctip = croot*taperRatio;
croot = plane.croot;
ctip = plane.ctip;

% find the chord distribution as a function of butt line
% chord = croot * (1 - buttLine) + ctip * buttLine;
% thickness = 0.9*tcRatio * chord;

% chord = croot * (1 - buttLine) + ctip * buttLine;

%% CENTER TUBE


centerTubeFraction = plane.centerTubeFraction;
plane.centerTubeSpan = plane.centerTubeFraction * b;

thicknessRoot = croot * tcRatio;
thicknessTip = ctip * tcRatio;

thicknessTransition = thicknessTip + (1 - centerTubeFraction) * ...
    (thicknessRoot - thicknessTip) / (1 - plane.centerSectionFraction);

% chordTransition = croot * (1 - centerTubeSpanFraction) + ctip * centerTubeSpanFraction
% thicknessTransition = tcRatio * chordTransition

tubeSkinClearance = 0.005; % spacing/clearance between skin and tube for ribs
centerTubeOuterDiameter = thicknessTransition - 2*tubeSkinClearance;
centerTubeOuterRadius = centerTubeOuterDiameter/2;

massEachSide = GMloadKG / (2 * plane.numberOfTubes); % each wing section takes half the load on GM
loadEachSide = massEachSide * 9.81; % force each wingtip must support, N

shear = ones(1, numPoints) * loadEachSide; % N
moment = shear.*(semispan-x); % N*m
momentRoot = moment(1);

% fudgefactor = 0.7333; % testing-based correction % from original Ben Boggs SparSizer
fudgeFactor = 1;

% ult_cap_stress_psi = 80000 * fudgefactor; % psi; 80000 for carbon uni (use weaker of tensile & compressive)
% CF data from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
ult_cap_stress = 850e6*fudgeFactor; % weaker of tensile and compressive strength since tubes/spars equally loaded
FS = 1.2; % factor of safety, 1.5 is std for ADT
stressAllowable = ult_cap_stress / FS;
plyThickness = 0.008*0.0254; % single ply thickness, converting from in --> m

centerTubeInnerRadius = (centerTubeOuterRadius^4 - ...
    (4/pi)*momentRoot*centerTubeOuterRadius/stressAllowable)^0.25;

centerTubeThickness = centerTubeOuterRadius - centerTubeInnerRadius;

centerTubeUniLayers = ceil(centerTubeThickness/plyThickness);
biaxialLayers = 4;
centerTubeLayers = centerTubeUniLayers + biaxialLayers;

centerTubeInnerRadius = centerTubeOuterRadius - centerTubeLayers*plyThickness;
centerTubeThickness = centerTubeOuterRadius - centerTubeInnerRadius;
successiveRadiiCenter = centerTubeOuterRadius:-plyThickness:centerTubeInnerRadius;
successiveCircumferencesCenter = successiveRadiiCenter * 2 * pi;
centerTubePlyArea = sum(successiveCircumferencesCenter) * centerTubeFraction * b;


if ~isreal(centerTubeInnerRadius) || centerTubeInnerRadius < 0
    fprintf('Center tube radius is imaginary or negative\n')
    plane.validTubes = false;
else
    plane.centerTubeOuterDiameter = centerTubeOuterDiameter;
    plane.centerTubeInnerDiameter = 2*centerTubeInnerRadius;
    plane.centerTubeThickness = centerTubeThickness;
    plane.centerTubeLayers = centerTubeLayers;
end



%% OUTER TUBE

tubeSemispanOverlap = 0.2;
outerTubeSpanFraction = (1 - centerTubeFraction + tubeSemispanOverlap)/2; % point of span where outer tube is loaded
outerTubeSpan = outerTubeSpanFraction * b;
tubeTubeGap = 0.002; % spacing between inner wall of center tube and outer wall of outboard tube

momentOuterTubeRoot = loadEachSide * outerTubeSpan;

if thicknessTip/2 - tubeSkinClearance < centerTubeInnerRadius - tubeTubeGap
    fprintf('Tip thickness limiting outer tube diameter\n')
end

outerTubeOuterRadius = min([centerTubeInnerRadius - tubeTubeGap, thicknessTip/2 - tubeSkinClearance]);
outerTubeRootInnerRadius = (outerTubeOuterRadius^4 - ...
    (4/pi)*momentOuterTubeRoot*outerTubeOuterRadius/stressAllowable)^0.25; % constant due to foam mandrel


outerPointXIndices = find(x > (1-2*outerTubeSpanFraction)*plane.b);
outerPointX = x(outerPointXIndices);
outerMoments = moment(outerPointXIndices);

% taking assumption of constant OUTER radius outer tubes
outerTubeInnerRadius = (outerTubeOuterRadius^4 - ...
    (4/pi)*outerMoments*outerTubeOuterRadius/stressAllowable).^0.25;

outerTubeThickness = outerTubeOuterRadius - outerTubeInnerRadius;

outerTubeUniLayers = ceil(outerTubeThickness/plyThickness);
outerTubeLayers = outerTubeUniLayers + biaxialLayers;

outerTubeInnerRadius = outerTubeOuterRadius - outerTubeLayers*plyThickness;

outerTubeThickness = outerTubeOuterRadius - outerTubeInnerRadius;

outerTubeAnnulusAreas = pi*outerTubeOuterRadius.^2 - pi*outerTubeInnerRadius.^2;

outerTubePlyVolume = sum( outerTubeAnnulusAreas(1:end-1) .* diff(outerPointX) );

outerTubePlyArea = outerTubePlyVolume/plyThickness;

tubePlyArea = plane.numberOfTubes * (centerTubePlyArea + 2*outerTubePlyArea);

%}





%{

b = plane.b
tcRatio = plane.tcRatio
GMloadKG = plane.GMloadingMultiplier * plane.seed_m2 % mass, in kg, put on plane during GM

% inputs
semispan = b / 2; % [in] half-span (single panel)

ds = 0.01; % numerical increment, m
x = (0:ds:semispan);
numPoints = length(x);
% buttLine = x/semispan; % fraction of way between root and tip

% cavg = 0.5*(croot + ctip)
% taperRatio = ctip/croot
% ctip = croot*taperRatio
% cavg = 0.5*(croot + croot*taperRatio)
% cavg = 0.5*croot*(1+taperRatio)
% croot = 2*cavg / (1 + taperRatio)
% croot = 2*cavg / (1 + taperRatio);
% ctip = croot*taperRatio;
croot = plane.croot
ctip = plane.ctip

% find the chord distribution as a function of butt line
% chord = croot * (1 - buttLine) + ctip * buttLine;
% thickness = 0.9*tcRatio * chord;

% chord = croot * (1 - buttLine) + ctip * buttLine;

%% CENTER TUBE


centerTubeFraction = plane.centerTubeFraction;
plane.centerTubeSpan = plane.centerTubeFraction * b;

thicknessRoot = croot * tcRatio
thicknessTip = ctip * tcRatio

thicknessTransition = thicknessTip + (1 - centerTubeFraction) * ...
    (thicknessRoot - thicknessTip) / (1 - plane.centerSectionFraction)

% chordTransition = croot * (1 - centerTubeSpanFraction) + ctip * centerTubeSpanFraction
% thicknessTransition = tcRatio * chordTransition

tubeSkinClearance = 0.005; % spacing/clearance between skin and tube for ribs
centerTubeOuterDiameter = thicknessTransition - 2*tubeSkinClearance
centerTubeOuterRadius = centerTubeOuterDiameter/2

massEachSide = GMloadKG / (2 * plane.numberOfTubes) % each wing section takes half the load on GM
loadEachSide = massEachSide * 9.81 % force each wingtip must support, N

shear = ones(1, numPoints) * loadEachSide; % N
moment = shear.*(semispan-x); % N*m
momentRoot = moment(1)

% fudgefactor = 0.7333; % testing-based correction % from original Ben Boggs SparSizer
fudgeFactor = 1;

% ult_cap_stress_psi = 80000 * fudgefactor; % psi; 80000 for carbon uni (use weaker of tensile & compressive)
% CF data from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
ult_cap_stress = 850e6*fudgeFactor; % weaker of tensile and compressive strength since tubes/spars equally loaded
FS = 1.5; % factor of safety, 1.5 is std for ADT
stressAllowable = ult_cap_stress / FS
plyThickness = 0.008*0.0254 % single ply thickness, converting from in --> m

centerTubeInnerRadius = (centerTubeOuterRadius^4 - ...
    (4/pi)*momentRoot*centerTubeOuterRadius/stressAllowable)^0.25

centerTubeThickness = centerTubeOuterRadius - centerTubeInnerRadius

centerTubeUniLayers = ceil(centerTubeThickness/plyThickness)
biaxialLayers = 4;
centerTubeLayers = centerTubeUniLayers + biaxialLayers

centerTubeInnerRadius = centerTubeOuterRadius - centerTubeLayers*plyThickness
centerTubeThickness = centerTubeOuterRadius - centerTubeInnerRadius
successiveRadiiCenter = centerTubeOuterRadius:-plyThickness:centerTubeInnerRadius
successiveCircumferencesCenter = successiveRadiiCenter * 2 * pi
centerTubePlyArea = sum(successiveCircumferencesCenter) * centerTubeFraction * b


if ~isreal(centerTubeInnerRadius) || centerTubeInnerRadius < 0
    fprintf('Center tube radius is imaginary or negative\n')
    plane.validTubes = false;
else
    plane.centerTubeOuterDiameter = centerTubeOuterDiameter;
    plane.centerTubeInnerDiameter = 2*centerTubeInnerRadius;
    plane.centerTubeThickness = centerTubeThickness;
    plane.centerTubeLayers = centerTubeLayers;
end



%% OUTER TUBE

outerTubeSpanFraction = (1 - centerTubeFraction + 0.15)/2 % point of span where outer tube is loaded
outerTubeSpan = outerTubeSpanFraction * b
tubeTubeGap = 0.002 % spacing between inner wall of center tube and outer wall of outboard tube

momentOuterTubeRoot = loadEachSide * outerTubeSpan;

if thicknessTip/2 - tubeSkinClearance < centerTubeInnerRadius - tubeTubeGap
    fprintf('Tip thickness limiting outer tube diameter\n')
end

outerTubeOuterRadius = min([centerTubeInnerRadius - tubeTubeGap, thicknessTip/2 - tubeSkinClearance])
outerTubeRootInnerRadius = (outerTubeOuterRadius^4 - ...
    (4/pi)*momentOuterTubeRoot*outerTubeOuterRadius/stressAllowable)^0.25 % constant due to foam mandrel


outerPointXIndices = find(x > (1-2*outerTubeSpanFraction)*plane.b)
outerPointX = x(outerPointXIndices);
outerMoments = moment(outerPointXIndices)

% taking assumption of constant OUTER radius outer tubes
outerTubeInnerRadius = (outerTubeOuterRadius^4 - ...
    (4/pi)*outerMoments*outerTubeOuterRadius/stressAllowable).^0.25

outerTubeThickness = outerTubeOuterRadius - outerTubeInnerRadius

outerTubeUniLayers = ceil(outerTubeThickness/plyThickness)
outerTubeLayers = outerTubeUniLayers + biaxialLayers

outerTubeInnerRadius = outerTubeOuterRadius - outerTubeLayers*plyThickness

outerTubeThickness = outerTubeOuterRadius - outerTubeInnerRadius

outerTubeAnnulusAreas = pi*outerTubeOuterRadius.^2 - pi*outerTubeInnerRadius.^2

outerTubePlyVolume = sum( outerTubeAnnulusAreas(1:end-1) .* diff(outerPointX) )

outerTubePlyArea = outerTubePlyVolume/plyThickness

tubePlyArea = plane.numberOfTubes * (centerTubePlyArea + 2*outerTubePlyArea)

%}

if ~isreal(outerTubeRootInnerRadius) || sum(outerTubeThickness < 0) > 0
    fprintf('Outer tube inner radius is imaginary or negative \n')
    plane.validTubes = false;
else
    plane.outerTubeOuterDiameter = 2*outerTubeOuterRadius;
    plane.outerTubeInnerDiameter = 2*outerTubeRootInnerRadius;
    plane.outerTubeLocations = outerPointX;
    plane.outerTubeThickness = outerTubeOuterRadius - outerTubeInnerRadius;
    plane.outerTubeLayers = outerTubeLayers;
    plane.outerTubeSpan = outerTubeSpan;
end


%% SHEAR WEB (FOAM)

shearWebVolume = plane.numberOfTubes * ( pi * centerTubeInnerRadius^2 * centerTubeFraction * b + ...
    sum(pi.*outerTubeInnerRadius(1:end-1).^2 .* diff(outerPointX)) );

if ~isreal(shearWebVolume)
    plane.validTubes = false;
end

tubeSizerSpreadsheetWeightFactor = 1.25;

tubePlyArea = tubePlyArea * tubeSizerSpreadsheetWeightFactor;
shearWebVolume = shearWebVolume * tubeSizerSpreadsheetWeightFactor;


%% CHECKS

if plane.validTubes == false % tubes are numerically impossible
    tubePlyArea = 2; % set ply area to high value (~= 1.5 kg) so weight is high and TOFL impossible
    shearWebVolume = 0.02;
    fprintf('Invalid Tubes')
end

end