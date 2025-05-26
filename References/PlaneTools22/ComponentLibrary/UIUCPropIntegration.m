clc;clear;close all;

path = "/Users/jailoonker/Documents/PlaneTools22/ComponentLibrary/UIUC-propDB/volume-4/data";
allFiles = dir(strcat(path));
allFiles(1:2,:) = [];

propDiams = [];
propPitchs = [];
dynamicDataFiles = [];

for i = 1:length(allFiles)
    numsInFile = regexp(allFiles(i).name,'\d*','Match');
    if length(numsInFile) == 4
        dynamicDataFiles = [dynamicDataFiles convertCharsToStrings(allFiles(i).name)];
        propDiams = [propDiams str2double(numsInFile(1))];
        propPitchs = [propPitchs str2double(numsInFile(2))];
    end
    propDiams = unique(propDiams);
    propPitchs = unique(propPitchs);
end

DVals = [];
PVals = [];
AVals = [];
BVals = [];
CVals = [];
EVals = [];
FVals = [];
GVals = [];
propNames = [];
for j = 1:length(propDiams)
    for k = 1:length(propPitchs)
        propellers = [];
        JVals = [];
        CTVals = [];
        CPVals = [];
        for i = 1:length(dynamicDataFiles) 
            numsInFile = regexp(dynamicDataFiles(i),'\d*','Match');
            if str2double(numsInFile(1)) == propDiams(j)
                if str2double(numsInFile(2)) == propPitchs(k)
                    % extract data here
                    % dataTable = dlmread(strcat(path, '/', dynamicDataFiles(i)), ' ', 1, 0)
                    dataTable = readmatrix(strcat(path, '/', dynamicDataFiles(i)));
                    JVals = [JVals; dataTable(:,1)];
                    CTVals = [CTVals; dataTable(:,2)];
                    CPVals = [CPVals; dataTable(:,3)];
                    
                end
            end
        end
        
        % curve fit and variable assignment code below
        if ~isempty(CTVals)
            CTfit = fit(JVals,CTVals,'poly2');
            ABC = coeffvalues(CTfit);
            CPfit = fit(JVals,CPVals,'poly2');
            EFG = coeffvalues(CPfit);
            
            DVals = [DVals; propDiams(j)];
            PVals = [PVals; propPitchs(k)];
            AVals = [AVals; ABC(1)];
            BVals = [BVals; ABC(2)];
            CVals = [CVals; ABC(3)];
            EVals = [EVals; EFG(1)];
            FVals = [FVals; EFG(2)];
            GVals = [GVals; EFG(3)];
            
        end
    end
end
 
allProps = [DVals PVals AVals BVals CVals EVals FVals GVals];
save('myFile2.mat','allProps');
