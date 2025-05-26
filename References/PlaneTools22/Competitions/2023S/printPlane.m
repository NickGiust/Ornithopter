function [] = printPlane(plane)

fprintf("\n--SCORES--\n")
    fprintf("Total Score: %.2f\n", plane.totalscore);
    fprintf("   M2 Score: %.2f\n", plane.m2results.m2score_normal); 
    fprintf("   M3 Score: %.2f\n", plane.m3results.m3score_normal); 

fprintf("\n--TAKE-OFF SPEEDS AND CLs--\n")
    fprintf("M2 VTO: %.2f\n", plane.m2results.VTO); 
    fprintf("M3 VTO: %.2f\n", plane.m3results.VTO);
    fprintf("---\n"); 
    fprintf("M2 CL_TO: %.2f\n", (plane.m_empty*9.81/plane.S) / (0.5*1.225*plane.m2results.VTO^2)); 
    fprintf("M3 CL_TO: %.2f\n", (plane.m_empty*9.81/plane.S) / (0.5*1.225*plane.m3results.VTO^2)); 

    fprintf("\n--M2 TAKEOFF W/ %.f%% THROTTLE (MAX 25 FT)--\n", 100*plane.m2results.minTOthrottle)
    fprintf("M2 TOFL: %s\n", plane.m2results.TOFLbool)
    if plane.m2results.TOFLbool == 'PASSED'
        fprintf("         %.2f ft\n", plane.m2results.TOFL/0.3048);
        TOsegmentIndex = find(plane.m2results.performance.V > plane.m2results.VTO, 1);
        fprintf("M2 Avg. Current, TO: %.2f A for %.2f s\n", mean(plane.m2results.performance.I(1:TOsegmentIndex-1)), plane.m2results.coursePoints.startTimes(2));
    end
    fprintf("\n--M3 TAKEOFF W/ %.f%% THROTTLE (MAX 25 FT)--\n", 100*plane.m3results.minTOthrottle)
    fprintf("M3 TOFL: %s\n", plane.m3results.TOFLbool)
    if plane.m3results.TOFLbool == 'PASSED'
        fprintf("         %.2f ft\n", plane.m3results.TOFL/0.3048);
        TOsegmentIndex = find(plane.m3results.performance.V > plane.m3results.VTO, 1);
        fprintf("M3 Avg. Current, TO: %.2f A for %.2f s\n\n", mean(plane.m3results.performance.I(1:TOsegmentIndex-1)), plane.m3results.coursePoints.startTimes(2));
    end

    if plane.m2results.TOFLbool == 'PASSED' 
        fprintf("M2 Max Prop Torque: %.2f Nm\n", max(plane.m2results.performance.proptorque));
    end
    if plane.m3results.TOFLbool == 'PASSED'
        fprintf("M3 Max Prop Torque: %.2f Nm\n", max(plane.m3results.performance.proptorque));
    end

    if plane.m2results.TOFLbool == "PASSED"
    fprintf("\n--LAPS/PERFORMANCE--\n")
    fprintf("M2: %1.f laps in %i minutes %.1f seconds\n    %.1f%% battery energy used\n    %.f%% throttle (optimized)\n", plane.m2results.nLaps, floor(plane.m2results.missionTime/60),rem(plane.m2results.missionTime,60),100*plane.m2results.batteryUsage,100*plane.m2results.optimizedThrottleSetting);
    fprintf("M2 Max. Cruise Speed: %.0f mph\n", max(plane.m2results.performance.V)*2.237)
    fprintf("M2 Cruise CL: %.2f\n", (plane.m_empty*9.81/plane.S) / (0.5*1.225*max(plane.m2results.performance.V)^2))
    end
    if plane.m3results.TOFLbool == "PASSED"
    fprintf("M3: %1.f laps in %.1f seconds\n    %.1f%% battery energy used\n", plane.m3results.nLaps, plane.m3results.missionTime, 100-100*plane.m3results.performance.s(end));
    fprintf("M3 Max. Cruise Speed: %.0f mph\n", max(plane.m3results.performance.V)*2.237)
    fprintf("M3 Cruise CL: %.2f\n", (plane.m_empty*9.81/plane.S) / (0.5*1.225*max(plane.m3results.performance.V)^2))
    end

    fprintf("\n--WING--\n")
    fprintf("Span: %.2f m\n", plane.b)
    fprintf("Chord: %.2f m\n", plane.c)
    fprintf("AR: %.1f\n", plane.AR)

    fprintf("\n--MASS--\n")
    fprintf("Empty Mass: %.2f kg\n", plane.m_empty);
    fprintf("M2 Mass: %.2f kg\n", plane.m2);
    fprintf("M3 Mass: %.2f kg\n", plane.m3);
    fprintf("Wing Mass: %.2f kg\n", plane.mWing);
    
    fprintf("\n--TAIL--\n")
    fprintf("Span: %.2f m\n", plane.bTail)
    fprintf("Chord: %.2f m\n", plane.cTail)
    fprintf("Height: %.2f m\n", plane.hTail)

    %% Print and Plot Mass Buildup
    % Initialize plotting arrays
    massPie = [plane.mFuse plane.mWing plane.mTail plane.mBat plane.mESC plane.mLg plane.mWiring plane.mMotor plane.mServos, plane.mPayloadMech, plane.mAvionics];
    massLabels = {'Fuselage', 'Wing', 'Tail', 'Battery', 'ESC', 'Landing Gear', 'Wiring', 'Motor', 'Servos', 'Payload Mech', 'Avionics'};
    
    % Print with loop through values and labels initialized above
    fprintf("\n--MASS BUILDUP--\n")
    for i=1:length(massPie)
        fprintf("%s: %.2f kg\n", massLabels{i}, massPie(i))
    end
    
    % Plot pie chart
    figure()
    ax = gca();
    p = pie(ax, massPie);
    ax.Colormap = turbo();
    pText = findobj(p,'Type','text');
    percentValues = get(pText,'String')';
    for i=1:length(pText)
        pText(i).String = sprintf("%s: %s\n(%.2f kg)", massLabels{i}, percentValues{i}, massPie(i));
    end

    legend(massLabels)

end