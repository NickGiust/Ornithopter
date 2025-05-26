function CheckLimits(mission, plane, missionResults)
performance = missionResults.performance;

if mission.missionNo == 2
    batMaxCurrent = plane.bat2maxCurrent*plane.nParallel2;
    D = plane.D2;
    fprintf('Mission 2 warnings:\n');
else
    batMaxCurrent = plane.bat3maxCurrent*plane.nParallel3;
    D = plane.D3;
    fprintf('\nMission 3 warnings:\n');
end

% General η for propeller stall
propEtaMin = .2;

%Exceeding battery current limit warning - ensures that I < Imax
if (max(performance.I) > batMaxCurrent)
    current = max(I)
    batMaxCurrent
    fprintf('Battery current warning: battery maximum discharge current exceeded.\n');
end

%Exceeding ESC current limit warning - ensures that I < Imax
if (max(performance.I) > plane.ESCMaxCurrent)
    current = max(performance.I)
    plane.ESCMaxCurrent
    fprintf('ESC current warning: ESC maximum current exceeded.\n');
end

%Exceeding motor power limit warning - ensures that Pelec < PelecMax
if (max(performance.Pelectric) > plane.motorMaxPower*plane.nMotors)
    power = max(performance.Pelectric)
    plane.motorMaxPower*plane.nMotors
    fprintf('Motor power warning: motor maximum power exceeded.\n');
end

%Prop stall warning - ensures that eta > etaMin
if max(performance.etaProp) < propEtaMin
    propEta = max(performance.etaProp)
    propEtaMin
    fprintf('Prop efficiency warning: propeller operating at very low efficiency and propeller stall likely.\n');
end

% Exceeding turn structural limit warning - ensures that n < nStruct
if max(performance.loadFactor) > plane.nStruct
    nTurn = max(performance.loadFactor)
    plane.nStruct
    fprintf('Structure warning: aircraft structural load factor limit exceeded during turn.\n');
end

% Cruise wing stall warning - ensures that CL < CLmax
if max(performance.CL) > plane.CLmax
    cruiseCL = max(performance.CL)
    plane.CLmax
    fprintf('Wing stall warning: CL exceeds wing maximum 3D CL.\n')
end

% Prop speed warning
propMachLimit = .9;
maxPropSpeed = propMachLimit*340*2/D;
if max(performance.propspeed) > maxPropSpeed
    propSpeed = max(performance.propspeed)
    maxPropSpeed
    fprintf('Propeller speed warning: prop blade tips exceed M = 0.7.\n')
end

end