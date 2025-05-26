function [area_sparCaps, shearWebVolume, plane] = sparSizerGM(plane, plotResults)

% Jack Ahrens - jahrens@usc.edu - 10/16/22
% based on sparSizer by Ben Boggs - boggsb@usc.edu - 10/30/18

if nargin < 2
    plotResults = false;
end

b = plane.b;
tcRatio = plane.tcRatio;
GMloadKG = plane.GMloadingMultiplier * plane.seed_m2; % mass, in kg, put on plane during GM

% inputs
semispan = b / 2; % [in] half-span (single panel)

ds = 0.01; % numerical increment, m
x = (0:ds:semispan);
numPoints = length(x);
buttLine = x/semispan; % fraction of way between root and tip

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
chord = croot * (1 - buttLine) + ctip * buttLine;
thickness = 0.9*tcRatio * chord;


massEachSide = GMloadKG / 2; % each wing section takes half the load on GM
loadEachSide = massEachSide * 9.81; % force each wingtip must support, N

shear = ones(1, numPoints) * loadEachSide; % N
moment = shear.*(semispan-x); % N*m

% bank_angle = acos(1/loadfactor)*180/pi; % unsuppress output to see max allowable bank angle

fudgefactor = 0.7333; % testing-based correction 

% use this link for carbon data! http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
ult_cap_stress_psi = 80000 * fudgefactor; % psi; 80000 for carbon uni (use weaker of tensile & compressive)
ult_cap_stress = ult_cap_stress_psi * 6894.76; % convert to Pa
FS = 1.3; % factor of safety, 1.5 is std for ADT
stressAllowable = ult_cap_stress / FS;
plyThickness = 0.008*0.0254; % single ply thickness, converting from in --> m


% define cap width as fraction of chord
% as measured in lab with calipers for UNCURED carbon
% should remeasure by laying up xx plies and dividing cured thickness by xx


% define cap width as fraction of chord - seems reasonable to use a cap
% with width 0.35*c for an airfoil like the e560
capChordFraction = 0.3;
capWidth = chord*capChordFraction; % user defined web thickness function. Spar cap will mimic

% assume inertia of spar caps is only due to the parallel axis theorem
% i.e. where the inertia of ONE rectangular cap section is usually:
% I = capWidth*capThickness^3/12 + capWidth*capThickness*(thickness/2)^2
% taking the summer plane 2023 case where:
% capWidth = 2; plyt = 0.0075; tcRatio = 0.14; chord = 0.21;
% the inertia I only decreases 4% by eliminating the first term, such that:
% I ~= capWidth*capThickness*(thickness/2)^2 for one of the two caps

capThicknessEach = ( moment.*(thickness/2)./ stressAllowable ) ./ (capWidth .* (thickness/2).^2 * 2);
numberOfPliesEach = ceil(capThicknessEach./plyThickness);
capThicknessAsBuilt = numberOfPliesEach.*plyThickness;

% calculate moment of inertia "as built"
I = capWidth.*capThicknessAsBuilt.*(thickness/2).^2 * 2;
failureMoment = ult_cap_stress*I./(thickness/2);

% multiply by two to account for both sides
area_sparCaps = 2 * trapz(x, numberOfPliesEach.*capWidth);
plane.wingSparCapPlies = numberOfPliesEach; % one side

% foam properties: https://www.dupont.com/content/dam/dupont/amer/us/en/...
% performance-building-solutions/public/documents/en/styrofoam-brand-...
% panel-core-30-pis-43-D100113-enUS.pdf

tau_foam = 35 * 6894.76; % convert 35 psi to Pa, = 241 kPa

% birch aircraft plywood properties:
% https://www.semanticscholar.org/paper/Comparison-of-Birch-and-Beech-...
% Wood-in-Terms-of-and-Cakiroglu-Demir/6d8b48cc36c3f8e954d26110ae918cdf980c9b0d

tau_aircraftPlywood = 3e6; % ~3 MPa, shear strength tau


% JAI: THIS SHOULD BE AN INPUT IN THE PLANE VARIABLE
webMaterial = plane.webMaterial;

if webMaterial == "plywood"
    tau = tau_aircraftPlywood;
elseif webMaterial == "foam"
    tau = tau_foam;
else
    error('incorrect shear web material entered')
end

% shear strength = shear force / web area

necessaryWebArea = loadEachSide / tau;

necessaryWebWidth = necessaryWebArea ./ thickness;

webWidth = max(necessaryWebWidth);

shearWebVolume = 2 * trapz(x, webWidth*thickness);

if plotResults
    
    %Bending Moment
    figure()
    hold on
    plot(x, moment, 'b', x, moment*FS, 'r--', x, failureMoment, 'k', 'LineWidth', 2);
    set(gca, 'FontName', 'Arial')
    ax = gca;
    ax.XLim = [0, semispan];
    %ax.YLim = [0, 7000];
    xticks(0:0.1:semispan);
    box off;
    legend('"Tip load" Moment', 'Allowable Moment (with FS)', 'Cap Failure Moment'); 
    set(gca,'FontSize',12);
    ylabel('Moment [N*m]');
    xlabel('Span Station [m]');


    disp('help')
    figure()
    subplot(4,1,1)
    plot(x, numberOfPliesEach, 'LineWidth', 2)
    ylabel('# plies per cap')
    set(gca, 'XTickLabel', [])


    disp('help 2')
    subplot(4,1,2)
    plot(x, capWidth*100, 'LineWidth', 2)
    ylabel('Cap width [cm]')
    ylim([0 1.1*max(capWidth*100)])
    set(gca, 'XTickLabel', [])

    subplot(4,1,3)
    plot(x, I*10^8, 'LineWidth', 2)
    ylabel('Cap inertia [cm^4]')
    

    subplot(4,1,4)
    plot(x, webWidth*ones(1, length(x))*100, 'LineWidth', 2)
    xlabel('Span Station [m]')
    ylabel('Web width [cm]')

    for i = 1:4
        subplot(4,1,i)
        xlim([0 semispan])
        set(gca, 'box', 'off')
        set(gca, 'FontName', 'Arial')
        set(gca,'FontSize',12);
    end


end

%{
% capw = Vcap*.0492*2; % old version, where Vcap was in in^3
% webw = Vweb*.0061*2; % old version, where Vcap was in in^3
sparw = capw + webw;
fprintf('Cap Weight (both sides): ' + string(capw) + ' lb (' + string(capw*16) + ' oz)\n');
fprintf('Web Weight (both sides): ' + string(webw) + ' lb (' + string(webw*16) + ' oz)\n');
fprintf('Spar Weight (both sides): ' + string(sparw) + ' lb (' + string(sparw*16) + ' oz)\n');
fprintf('Spar Weight (one side): ' + string(sparw/2) + ' lb (' + string(sparw*8) + ' oz)\n'); 
fprintf('\n')

%% Pieces to cut of uni-directional CF
fprintf('----------------------------------------------\n')
fprintf('PIECES TO CUT OF UNI-DIRECTIONAL CARBON FIBER:\n')
fprintf('\n')

fprintf('Spar cap width: ' + string(capWebWidth(1)) + ' in.\n')
ply_length_v = size(ply_length);
for ii = 1:ply_length_v(2)
    fprintf('Cut 2 plies of ' + string(ply_length(ii)) + ' in.\n'); 
end
fprintf('Total plies: ' + string(2*ply_length_v(2)) + '\n')

end
%}
