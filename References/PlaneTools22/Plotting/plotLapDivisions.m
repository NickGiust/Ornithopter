missionPoints = missionResults.missionPoints;
coursePoints = missionResults.coursePoints;
yDivide = .85;
yLapNo = .9;
ySegment = .5;
m3segmentLabels = 0;

colors = [255 0 0; ...
    255 127 0; ...
    255 255 0; ...
    0 255 0; ...
    0 0 255; ...
    75 0 130; ...
    148 0 211]./255;
colors = [colors; colors; colors; colors; colors];

hold on;
% Plot lap divisions
for lapNo = 1:length(missionPoints.lapNos)
    startTime = missionPoints.lapStartTimes(lapNo);
    try endTime = missionPoints.lapStartTimes(lapNo+1);
    catch endTime = max(xlim); end
    top = max(ylim);
    lapColor = [colors(lapNo,1) colors(lapNo,2) colors(lapNo,3)];
    bottom = min(ylim) + yDivide*(max(ylim)-min(ylim));
    
    patch([startTime endTime endTime startTime], [top top bottom bottom], ...
       lapColor, 'FaceAlpha', 0.3, 'LineStyle', ':');
    
    yText = min(ylim) + yLapNo*(max(ylim)-min(ylim));
    xText = startTime + .5*(endTime-startTime);
    text(xText,yText,num2str(lapNo),'FontUnits','pixels','FontSize',12,'HorizontalAlignment','center')
end
% Plot lap segment divisions
hAx = gca;
for segmentNo = 1:length(coursePoints.segments)
    startTime = coursePoints.startTimes(segmentNo);
    try endTime = coursePoints.startTimes(segmentNo+1);
    catch endTime = max(xlim); end
    top = min(ylim) + yDivide*(max(ylim)-min(ylim));
    bottom = min(ylim);
    
    patch([startTime endTime endTime startTime], [top top bottom bottom], ...
        [0 0 0], 'FaceAlpha', 0, 'LineStyle', ':');
    
    if mission.missionNo ~= 3 | m3segmentLabels
        yText = min(ylim) + ySegment*(max(ylim)-min(ylim));
        xText = startTime + .5*(endTime-startTime);
        t = text(xText,yText,coursePoints.segments(segmentNo),'FontUnits','normalized','FontSize',.05,'HorizontalAlignment','center');
        set(t,'Rotation',90);
    end
end
hold off;
%patch([0 1E+5 1E+5 0], [max(ylim) max(ylim) 0 0], 'r')
