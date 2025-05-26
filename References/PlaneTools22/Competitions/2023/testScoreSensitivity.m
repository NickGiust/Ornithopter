% for testing scoreSensitivity.m

clear; close all; clc;

%% Best plane for each mission

% M2 assumptions
B.mPackageBest = 7;
B.numLapsBest = 18;
B.M2Best = B.mPackageBest*B.numLapsBest;

% M3 assumptions
B.lengthBest = 1.1;
B.timeBest = 50; % for M3, fixed 3 laps
B.M3Best = B.lengthBest/B.timeBest;

B.GMBest = 50; % best team's GM loading multiplier


%% Our plane (USC)

X.numLapsUSC = 16; % fly this many laps with the package
X.mPackageUSC = 3.67; % m2*payloadFraction2

X.timeUSC = 67.6; % seconds
X.lengthUSC = 0.95; % meters

X.mTestUSC = 672;
X.mMaxUSC = 11.1; % max of M2 and M3 gross aircraft masses
% X.GMUSC = X.mTestUSC / X.mMaxUSC;
X.GMUSC = 79;

%% execute sensitivity analysis

% dSdX and dSdB variables represent the relative score sensitivity to each
% parameter

[dSdX_Base, dSdB_Base, S_base] = relativeScoreSensitivity(X, B)

changes = 0.6:0.00001:1.4;

XArray = X;
XfieldNames = fieldnames(X)

figure()
hold on
xlabel('Parameter Change [%]')
ylabel('Change in Mission Score [%]')

doNotPlotIndices = 7;

% legendList = strings(1,length(XfieldNames));

LineTypes = ["-", "--", "-.", ":"];

for i = 1:length(XfieldNames)
    currentField = XfieldNames{i};

    XArray.(currentField) = XArray.(currentField) * changes;

    if strcmpi(currentField, 'numLapsUSC')
        XArray.(currentField) = floor(XArray.(currentField));
    end

    [dSdX_Current, dSdB_Current, S_Current] = relativeScoreSensitivity(XArray, B, changes);
    
    relativeScoreChange = (S_Current - S_base) / S_base;

    if i ~= doNotPlotIndices

        LineTypeCurrent = LineTypes(mod(i-1,length(LineTypes))+1);
        plot(changes*100-100, relativeScoreChange*100, LineTypeCurrent, 'LineWidth', 2)

    end

    % legendList(i) = string(currentField(1:end-3));

    XArray.(currentField) = X.(currentField);
    %% INTERSECTIONS ON Y AXIS SHOULD BE ZERO
end

% legendList(doNotPlotIndices) = []

legendList = ["Mission 2 Laps", "Electronics Package Mass", "Mission 3 Time",...
    "Antenna Length", "Ground Mission Test Mass", "Maximum Aircraft Mass"];

plotCurrentDesignPoint = false;

if plotCurrentDesignPoint
    plot(0, 0, 'pentagram', 'MarkerSize', 25, 'MarkerFaceColor', 'k', 'LineWidth', 1, 'Color','k')
    % plot([-50 50], [0 0], 'k--', 'LineWidth', 2)
    legend([legendList, "Current Design Point"], 'Location', 'best')
else
    legend(legendList, 'Location', 'best')
end

% grid minor

grid on
ax = gca;
ax.FontSize = 14;