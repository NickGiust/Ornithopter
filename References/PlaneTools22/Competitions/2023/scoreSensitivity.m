function [dSdX, dSdB] = scoreSensitivity(X, B)
%%SCORESENSITIVITY finds the differential change in score with respect to
%%change in the score variables, X, like antenna length or package mass
% X Should contain (in structure form): mpackageUSC, numLapsUSC,
% lengthUSC, timeUSC, mTestUSC, and mMaxUSC 
% B represents the best score variables, and should contain:
% mpackageBest, numLapsBest, lengthBest, timeBest, and GMMultiplierBest
% Note: mTestUSC is the amount of weight loaded on the plane and
% mMaxUSC is the greater of the M2 and M3 weights (probably M2)


% Score = S = M2 + M3 + GM

% M2 = mPackageUSC*numLapsUSC / (mPackageBest*numLapsBest)
% M2 = mPackageUSC*numLapsUSC / M2Best

% M3 = (lengthUSC/timeUSC) / (lengthBest/timeBest)
% M3 = (lengthUSC/timeUSC) / M3Best
% M3 = lengthUSC*timeBest / (timeUSC*lengthBest)

% GM = (mTestUSC/mMaxUSC) / (mTestBest/mMaxBest)
% GM = (mTestUSC/mMaxUSC) / (GMMultiplierBest)


% find dS/dX, where X is the vector of each of the individual contributors to score

dSdX.mPackageUSC = X.numLapsUSC/B.M2Best - X.mTestUSC / (X.mMaxUSC^2*B.GMBest);
dSdX.numLapsUSC = X.mPackageUSC/B.M2Best;
dSdB.M2Best = - X.mPackageUSC*X.numLapsUSC / B.M2Best^2;
dSdB.mPackageBest = - X.mPackageUSC*X.numLapsUSC / (B.mPackageBest^2*B.numLapsBest);
dSdB.numLapsBest = - X.mPackageUSC*X.numLapsUSC / (B.mPackageBest*B.numLapsBest^2);


dSdX.lengthUSC = 1 / (X.timeUSC*B.M3Best);
dSdX.timeUSC = - X.lengthUSC / (X.timeUSC^2*B.M3Best);
dSdB.M3Best = - X.lengthUSC / (X.timeUSC*B.M3Best^2);
dSdB.lengthBest = - X.lengthUSC*B.timeBest / (X.timeUSC*B.lengthBest^2);
dSdB.timeBest = - X.lengthUSC / (X.timeUSC*B.lengthBest);

dSdX.mTestUSC = 1 / (X.mMaxUSC*B.GMBest);
dSdX.mMaxUSC = - X.mTestUSC / (X.mMaxUSC^2*B.GMBest);
dSdB.GMBest = - X.mTestUSC / (X.mMaxUSC*B.GMBest^2);
dSdB.GMUSC = 1 / B.GMBest; % chance in score w.r.to change in our GMloadingMultiplier



end
