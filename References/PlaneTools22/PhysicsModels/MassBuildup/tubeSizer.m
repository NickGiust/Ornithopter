function [tubePlyArea, shearWebVolume, plane] = tubeSizer(plane)

% Jack Ahrens - jahrens@usc.edu - 10/16/22

plane.validTubes = true;

% %{

b = plane.b;
% AR = plane.AR;
tcRatio = plane.tcRatio;
GMloadKG = plane.GMloadingMultiplier * plane.seed_m2; % mass, in kg, put on plane during GM
GMloadKG = 650;

% inputs
semispan = b / 2; % [in] half-span (single panel)

croot = plane.croot;
ctip = plane.ctip;



%% INPUT PARAMETERS

massEachSide = GMloadKG / (2 * plane.numberOfTubes); % each wing section takes half the load on GM
loadEachSide = massEachSide * 9.81; % force each wingtip must support, N

momentRoot = loadEachSide * semispan;

% fudgefactor = 0.7333; % testing-based correction % from original Ben Boggs SparSizer
fudgeFactor = 1;

% ult_cap_stress_psi = 80000 * fudgefactor; % psi; 80000 for carbon uni (use weaker of tensile & compressive)
% CF data from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
ult_cap_stress = 850e6*fudgeFactor; % weaker of tensile and compressive strength since tubes/spars equally loaded
FS = 1.3; % factor of safety, 1.5 is std for ADT
stressAllowable = ult_cap_stress / FS;

plyThickness = 0.008*0.0254; % single ply thickness, converting from in --> m

tubeSkinClearance = 0.005; % spacing/clearance between skin and tube for ribs
tubeTubeGap = 0.002; % spacing between inner wall of center tube and outer wall of outboard tube


biaxialLayers = 4;


%% OUTER TUBE



outerSectionSpanFraction = plane.outerSectionSpan / semispan;
plane.outerSectionSpanFraction = outerSectionSpanFraction;

effectiveCenterSectionSemispan = semispan - plane.outerSectionSpan; % centerline to outer wing section

spansOuterSection = effectiveCenterSectionSemispan:0.01:semispan;

buttLineOuter = linspace(0, 1, length(spansOuterSection));

% find the chord distribution as a function of butt line
outerChords = croot * (1 - buttLineOuter) + ctip * buttLineOuter;
outerThicknesses = 0.95 * tcRatio * outerChords;

outerTubeDiameters = outerThicknesses - 2 * tubeSkinClearance;

outerTubeOuterRadii = outerTubeDiameters / 2;

numPointsOuter = length(spansOuterSection);

shearOuter = ones(1, numPointsOuter) * loadEachSide; % N
momentsOuter = shearOuter.*(semispan-spansOuterSection);

outerTubeInnerRadii = (outerTubeOuterRadii.^4 - ...
    (4/pi).*momentsOuter.*outerTubeOuterRadii/stressAllowable).^0.25;


outerTubeThicknesses = outerTubeOuterRadii - outerTubeInnerRadii;

outerTubeUniLayers = ceil(outerTubeThicknesses/plyThickness);
outerTubeLayers = outerTubeUniLayers + biaxialLayers;

outerTubeInnerRadii = outerTubeOuterRadii - outerTubeLayers*plyThickness;

outerTubeThicknesses = outerTubeOuterRadii - outerTubeInnerRadii;

outerTubeAnnulusAreas = pi*outerTubeOuterRadii.^2 - pi*outerTubeInnerRadii.^2;

outerTubePlyVolume = sum( outerTubeAnnulusAreas(1:end-1) .* diff(spansOuterSection) );

outerTubePlyArea = outerTubePlyVolume/plyThickness;




%% CENTER TUBE

tubeOverlapFraction = 0.2; % fraction of semispan that center and outer tubes overlap
centerTubeSpan = 2 * (effectiveCenterSectionSemispan + tubeOverlapFraction * semispan);

outerTubeInnerRadius = min(outerTubeInnerRadii(spansOuterSection < spansOuterSection(1) + tubeOverlapFraction*semispan));

centerTubeOuterDiameter = 2 * ( outerTubeInnerRadius - tubeTubeGap );


% centerTubeOuterDiameter = 2 * ( outerTubeInnerRadius - tubeTubeGap );
centerTubeOuterRadius = centerTubeOuterDiameter/2;



centerTubeInnerRadius = (centerTubeOuterRadius^4 - ...
    (4/pi)*momentRoot*centerTubeOuterRadius/stressAllowable)^0.25;

centerTubeThickness = centerTubeOuterRadius - centerTubeInnerRadius;

centerTubeUniLayers = ceil(centerTubeThickness/plyThickness);
centerTubeLayers = centerTubeUniLayers + biaxialLayers;

centerTubeInnerRadius = centerTubeOuterRadius - centerTubeLayers*plyThickness;
centerTubeThickness = centerTubeOuterRadius - centerTubeInnerRadius;
successiveRadiiCenter = centerTubeOuterRadius:-plyThickness:centerTubeInnerRadius;
successiveCircumferencesCenter = successiveRadiiCenter * 2 * pi;
centerTubePlyArea = sum(successiveCircumferencesCenter) * centerTubeSpan;



%% Check tube validity, pass results out


if ~isreal(outerTubeInnerRadii) || sum(outerTubeThicknesses < 0) > 0 || outerTubeInnerRadius <= 0
    fprintf('Outer tube inner radius is imaginary or negative \n')
    plane.validTubes = false;
else
    plane.outerTubeInnerDiameter = 2*outerTubeInnerRadius;
    plane.outerTubeLocations = spansOuterSection;
    plane.outerTubeThickness = outerTubeOuterRadii - outerTubeInnerRadii;
    plane.outerTubeLayers = outerTubeLayers;
end


if ~isreal(centerTubeOuterRadius) || ~isreal(centerTubeInnerRadius) || centerTubeInnerRadius < 0 || centerTubeOuterDiameter <= 0
    fprintf('Center tube radius is imaginary or negative\n')
    plane.validTubes = false;
else
    plane.centerTubeSpan = centerTubeSpan;
    plane.centerTubeSpanFraction = centerTubeSpan / b;
    plane.centerTubeOuterDiameter = centerTubeOuterDiameter;
    plane.centerTubeInnerDiameter = 2*centerTubeInnerRadius;
    plane.centerTubeThickness = centerTubeThickness;
    plane.centerTubeLayers = centerTubeLayers;
end




%% CHECKS

if plane.validTubes == false % tubes are numerically impossible
    tubePlyArea = 2; % set ply area to high value (~= 1.5 kg) so weight is high and TOFL impossible
    shearWebVolume = 0.02;
    fprintf('Invalid Tubes')

else

    tubePlyArea = plane.numberOfTubes * (centerTubePlyArea + 2*outerTubePlyArea);

    shearWebVolume = plane.numberOfTubes * ( pi * centerTubeInnerRadius^2 * centerTubeSpan + ...
    sum(pi.*outerTubeInnerRadii(1:end-1).^2 .* diff(spansOuterSection)) );

    if ~isreal(shearWebVolume)
        plane.validTubes = false;
    end
    
    tubeSizerSpreadsheetWeightFactor = 1.2;
    
    tubePlyArea = tubePlyArea * tubeSizerSpreadsheetWeightFactor;
    shearWebVolume = shearWebVolume * tubeSizerSpreadsheetWeightFactor;

end


end