clear;
load tradestudydata.mat
noEggs = 1:21;
chord = linspace(0,12,40);

%contourf(noEggs,chord,datamat2)                    %uses contour lines
surf(noEggs,chord,datamat2,'FaceColor','interp')    %just shading, no contour lines

view(2)                 %places the view straight down (rather than isometric)
axis([1 21 0 12])       %sets axis limits
shading interp          %blends the shading
colorbar                %adds the color scale for z
xlabel('No. of eggs')
ylabel('Chord [in]')
title('Score vs. no of eggs and chord')