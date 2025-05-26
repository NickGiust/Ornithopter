function [m2scores_normal,m3scores_normal] = normalizeScores(m2scores_USC,m3scores_USC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% M2

% Best plane assumptions
radarMass_best = 8;
nLaps = 15;
m2score_best = radarMass_best * nLaps;

% Normalize
m2scores_normal = m2scores_USC / m2score_best;
m2scores_normal(m2scores_normal > 1) = 1;

% M3
time_best = 45;
m3score_best = 3 / time_best;

m3scores_normal = m3scores_USC/m3score_best;
m3scores_normal(m3scores_normal > 1) = 1;

end

