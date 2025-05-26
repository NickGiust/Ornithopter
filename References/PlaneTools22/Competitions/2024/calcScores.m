function [m2scores,m3scores,totalscores,planes] = calcScores(planes)
%% Summary

% calcScores takes the results of the missions from SimulateMission and
% applies the year's scoring equations to output the official score of the
% plane. Best score assumptions are also made here

%% Best plane assumptions

% M2
massPayload2_best = 3; % [kg]
time2_best = 66; % [s]
m2score_best = massPayload2_best / time2_best;

% M3
nPassengers_best = 75; % [-]
nLaps_best = 12; % [-]
batCap_best = 2.2; % [Ah]
m3score_best = nPassengers_best * nLaps_best / batCap_best;

%% Score calculation

% preallocation
numPlanes = length(planes);
m2scores = zeros(1, numPlanes);
m3scores = zeros(1, numPlanes);

for i = 1:numPlanes
    % m2 scores
    if planes{i}.m2results.missionTime > 0 % only give a score if the time is greater than 0
        m2score_USC = planes{i}.massPayloads2 / planes{i}.m2results.missionTime; % competition rule
    else
        m2score_USC = 0; % set score = 0 if the time is <= 0 (meaning mission not run)
    end
    m2score_normal = m2score_USC / m2score_best; % normalization
    m2score_normal(m2score_normal > 1) = 1; % if better than best prediction, set score to 1
    if m2score_normal > 0
        m2score_total = m2score_normal + 1; % add completion constant to score
    else
        m2score_total = m2score_normal;
    end
    
    planes{i}.m2results.m2score_normal = m2score_total;
    m2scores(i) = m2score_total;

    % m3 scores
    if planes{i}.m3results.missionTime > 0 % only give a score if the time is greater than 0
        m3score_USC = planes{i}.m3results.nLaps * planes{i}.nPayloads3 / planes{i}.bat3capacity; % update to battery capacity
    else
        m3score_USC = 0; % set score = 0 if the time is <= 0 (meaning mission not run)
    end
    m3score_normal = m3score_USC / m3score_best;
    m3score_normal(m3score_normal > 1) = 1;
    if m3score_normal > 0
        m3score_total = m3score_normal + 2; % add completion constant to score
    else
        m3score_total = m2score_normal;
    end
    
    planes{i}.m3results.m3score_normal = m3score_total;
    m3scores(i) = m3score_total;

    planes{i}.totalscore = m2score_total + m3score_total;
end

totalscores = m2scores + m3scores;

end