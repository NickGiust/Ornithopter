% Ben Boggs - boggsb@usc.edu - 10/30/18

clear all; clc; close all;

% inputs
span = 48; % [in] half-span (single panel)
chord = 10; % [in]
tc = 0.14; % t/c
weight = 18; % max weight plane will fly with in [lb]
loadfactor = 4; % max g's 

bank_angle = acos(1/loadfactor)*180/pi; % unsuppress output to see max allowable bank angle

fudgefactor = .7333; % testing-based correction 

ult_cap_stress = 80000 * fudgefactor; % psi; 80000 for carbon uni (use weaker of tensile & compressive)
safetyf = 1.5; % factor of safety, 1.5 is std for ADT
plyt = .0075; % single ply thickness in [in].
ds = .01; % numerical increment

x = (0:ds:span);

capWebWidth = 0.*x + 2; % user defined web thickness function. Spar cap will mimic

avgpanellift = weight*loadfactor/2;
semiarea = span*chord;
shear = avgpanellift-(avgpanellift.*x/span);
moment = shear.*(span-x)./2;

% calculations assume balsa shear web 

Vcap = 0;
Vweb = 0;
prevplies = 0;
cplies = 1;
for (i=1:length(x))
    mall(i) = moment(i)*safetyf;
    I = 2 * ( (capWebWidth(i)/12)*plyt^3 + (capWebWidth(i)*plyt*( ( (tc*chord) - plyt)/2)^2) );
    M = ult_cap_stress*I/(tc*chord/2);
    j = 1;
    cplies = 1;
    while M < mall
        j = j + 1;
        cplies = cplies + 1;
        b = capWebWidth(i);
        h = plyt*j;
        d = (tc*chord-h)/2;
        I = 2 * ( ((b*h^3)/12) + b*h*d^2 );
        M = ult_cap_stress*I/(tc*chord/2);
    end
    
    if (cplies ~= prevplies) 
%         if (cplies == 1)
%             str = ' ply at x = ';
%         else
%             str = ' plies at x = ';
%         end
%         fprintf(string(cplies) + str + string(x(i)) + ' in.\n');
        %% Spar Caps Length Calculations (to print)
        ply_length(cplies) = 2*x(i);
        if cplies == 1
            ply_length = circshift(ply_length,1);
            ply_length(cplies) = 2 * span; % [in] (2 * half-span)
        end
    end
    prevplies = cplies;
    
    Vcap = Vcap + ds*j*plyt*capWebWidth(i)*2;
    Vweb = Vweb + (tc*chord*capWebWidth(i)*ds)-(ds*j*plyt*capWebWidth(i)*2);
    cap(i) = j;
    capm(i) = M;
end

%Bending Moment
figure(1)
hold on;
plot(x, moment, 'b', x, mall, 'r--', x, capm, 'k', 'LineWidth', 2);
set(gca, 'FontName', 'Arial')
ax = gca;
ax.XLim = [0, span];
%ax.YLim = [0, 7000];
xticks((0:6:span));
box off;
legend('Lifting Moment', 'Allowable Moment', 'Cap Failure Moment'); 
set(gca,'FontSize',12);
ylabel('Moment [lb*in]');
xlabel('Span Station [in]');

capw = Vcap*.0492*2;
webw = Vweb*.0061*2;
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