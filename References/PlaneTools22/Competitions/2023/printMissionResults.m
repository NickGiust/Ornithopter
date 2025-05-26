function printMissionResults(missionResults, mass)
%PRINTMISSIONRESULTS Summary of this function goes here
%   Detailed explanation goes here
fprintf('---------------- MISSION %i RESULTS ----------------\n', missionResults.missionNo)
fprintf('Mass: %.2f kg\n', mass)
fprintf('TOFL: %.1f ft\n', missionResults.TOFL)
fprintf('# of Laps: %i\n', missionResults.nLaps)
fprintf('Completion Time: %.1f sec\n', missionResults.missionTime)
fprintf('Score: %.2f\n', missionResults.score)
end

