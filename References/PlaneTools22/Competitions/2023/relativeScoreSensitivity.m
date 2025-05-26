function [dSdX, dSdB, S] = relativeScoreSensitivity(X, B, changes)
%%SCORESENSITIVITY finds the differential change in score with respect to
%%change in the score variables, X, like antenna length or package mass
% X Should contain (in structure form): mpackageUSC, numLapsUSC,
% lengthUSC, timeUSC, mTestUSC, and mMaxUSC 
% B represents the best score variables, and should contain:
% mpackageBest, numLapsBest, lengthBest, timeBest, and GMMultiplierBest
% Note: mTestUSC is the amount of weight loaded on the plane and
% mMaxUSC is the greater of the M2 and M3 weights (probably M2)

if nargin < 3
    changes = 1;
end

if sum(changes==1) == 0
    error('Need to have a ''base'' value where the parameter is unchanged')
end

if length(X.mPackageUSC) > 1

    mPackageUSC_base = X.mPackageUSC(changes == 1);
    
    X.mMaxUSC = X.mMaxUSC - mPackageUSC_base + X.mPackageUSC;

end

M2 = X.mPackageUSC.*X.numLapsUSC ./ (B.mPackageBest.*B.numLapsBest);

% M2 = mPackageUSC*numLapsUSC / (mPackageBest*numLapsBest)
% M2 = mPackageUSC*numLapsUSC / M2Best

M3 = (X.lengthUSC./X.timeUSC) ./ (B.lengthBest./B.timeBest);

% M3 = (lengthUSC/timeUSC) / (lengthBest/timeBest)
% M3 = (lengthUSC/timeUSC) / M3Best
% M3 = lengthUSC*timeBest / (timeUSC*lengthBest)

GM = (X.mTestUSC./X.mMaxUSC) ./ (B.GMBest);

% GM = (mTestUSC/mMaxUSC) / (mTestBest/mMaxBest)
% GM = (mTestUSC/mMaxUSC) / (GMMultiplierBest)

S = M2 + M3 + GM; % total score


% find (dS/dX) * (X/S), where X is the vector of each of the individual contributors to score

dSdX.mPackageUSC = ( X.numLapsUSC./B.M2Best - X.mTestUSC ./ (X.mMaxUSC.^2.*B.GMBest) ) .* (X.mPackageUSC./S);
dSdX.numLapsUSC = (X.mPackageUSC./B.M2Best) .* (X.numLapsUSC./S);
dSdB.M2Best = - (X.mPackageUSC.*X.numLapsUSC / B.M2Best^2) .* (B.M2Best./S);
dSdB.mPackageBest = - (X.mPackageUSC.*X.numLapsUSC / (B.mPackageBest.^2*B.numLapsBest)) .* (B.mPackageBest./S);
dSdB.numLapsBest = - (X.mPackageUSC.*X.numLapsUSC / (B.mPackageBest.*B.numLapsBest.^2)) .* (B.numLapsBest./S);


dSdX.lengthUSC = (1 ./ (X.timeUSC.*B.M3Best)) .* (X.lengthUSC./S);
dSdX.timeUSC = - (X.lengthUSC ./ (X.timeUSC.^2.*B.M3Best)) .* (X.timeUSC./S);
dSdB.M3Best = - (X.lengthUSC ./ (X.timeUSC.*B.M3Best.^2)) .* (B.M3Best./S);
dSdB.lengthBest = - (X.lengthUSC.*B.timeBest ./ (X.timeUSC.*B.lengthBest.^2)) .* (B.lengthBest./S);
dSdB.timeBest = - (X.lengthUSC ./ (X.timeUSC.*B.lengthBest)) .* (B.timeBest./S);

dSdX.mMaxUSC = - (X.mTestUSC ./ (X.mMaxUSC.^2.*B.GMBest)) .* (X.mMaxUSC./S);
dSdX.mTestUSC = (1 ./ (X.mMaxUSC.*B.GMBest)) .* (X.mTestUSC./S);
dSdX.GMUSC = (1 ./ B.GMBest) .* (X.GMUSC./S); % chance in score w.r.to change in our GMloadingMultiplier
dSdB.GMBest = - (X.mTestUSC ./ (X.mMaxUSC.*B.GMBest.^2)) .* (B.GMBest./S);



end
