function plotPerformance(plane, missionNo, yName, holdGraph, lineStyleCode, secondAxis, userLineWidth)

if nargin < 4
    holdGraph = false; % assume hold is off unless otherwise specified
end

if nargin < 5
    lineStyleCode = '-'; % assume solid line unless otherwise specified
end

if nargin < 6
    secondAxis = false;
end

if nargin < 7
    userLineWidth = 2;
end

% sort based on mission
if missionNo == 2
    results = 'm2results';
elseif missionNo == 3
    results = 'm3results';
else
    return
end

if ~isfield(plane.(results), "performance")
    fprintf('Couldn''t plot %s because mission %.0f not run\n', yName, missionNo)
    return
end

% check that y variable given is in the table
if ~isfield(plane.(results).performance, yName)
    fprintf('Given field for plotPerformance is not calculated or stored\n')
    return
end

% getting x and y values
t = plane.(results).performance.t;
y = plane.(results).performance.(yName);

% setting up colors
colorPalette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9"];
persistent colorIndex
if isempty(colorIndex)  % Initialize colorIndex if it's empty
    colorIndex = 1;
end

% setting up plot
if ~holdGraph
    figure()
    hold on
    colorIndex = 1; % reset color to first one if using a new graph
end

% adding labels
if ~holdGraph
    titleText = sprintf('Mission: %.0f', missionNo);
    title(titleText)
    xlabelText = sprintf('t [s]');
    xlabel(xlabelText)
    ylabelText = sprintf('%c', yName);
    ylabel(ylabelText)

    % separate segments (with labels)
%     for i = 3:length(plane.(results).coursePoints.startTimes)
%         if strcmpi(plane.(results).coursePoints.segments(i), 'halfStraightaway')
%             xline(plane.(results).coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Straightaway')
%         elseif strcmpi(plane.(results).coursePoints.segments(i), 'turn180')
%             xline(plane.(results).coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Turn')
%         elseif strcmpi(plane.(results).coursePoints.segments(i), 'turn360')
%             xline(plane.(results).coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', '360')
%         elseif strcmpi(plane.(results).coursePoints.segments(i), 'takeoff')
%             xline(plane.(results).coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Takeoff')
%         elseif strcmpi(plane.(results).coursePoints.segments(i), 'climb')
%             xline(plane.(results).coursePoints.startTimes(i), 'LineStyle', '--', 'Color', '#d3d3d3', 'Label', 'Climb')
%         end
%     end

elseif holdGraph && secondAxis
    yyaxis right
    ylabelRightText = sprintf('%c', yName);
    ylabel(ylabelRightText)
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
end

% plotting
plot(t, y, 'Color', colorPalette(colorIndex), 'LineWidth', userLineWidth, 'LineStyle', lineStyleCode)

% incrementing the color palette
colorIndex = colorIndex + 1;
if colorIndex > length(colorPalette)
    colorIndex = 1;
end

end