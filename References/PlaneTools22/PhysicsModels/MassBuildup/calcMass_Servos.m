function [massServo] = calcMass_Servos(servo_torque, servo_voltage)

% INPUT     servo_torque - torque of servo used
% OUTPUT    servo_mass - mass of servo

% Sample data - {Name, Voltage, torque [kg-cm], mass [g]}
servo_data_lowV = {'HS-485HB', 4.8, 4.82, 45.1; ...
                   'HS-425BB', 4.8, 3.31, 45.4; ...
                   'HS-85MG', 4.8, 3.02, 21.8; ...
                   'HS-325HB', 4.8, 3.02, 42.8; ...
                   'HS-45HB', 4.8, 1.01, 7.9; ...
                   'HS-525BB', 4.8, 3.30, 45.2; ...
                   'HS-635HB', 4.8, 4.97, 49.9; ...
                   'HS-805MG', 4.8, 19.8, 197; ...
                   'HS-81', 4.8, 2.59, 16.4; ...
                   'HS-945MG', 4.8, 8.80, 55.9};
               
servo_data_highV = {'HS-485HB', 6.0, 5.98, 45.1; ...
                    'HS-425BB', 6.0, 4.10, 45.4; ...
                    'HS-85MG', 6.0, 3.53, 21.8; ...
                    'HS-325BB', 6.0, 3.67, 42.8; ...
                    'HS-45HB', 6.0, 1.22, 7.9; ...
                    'HS-525BB', 6.0, 4.10, 45.2; ...
                    'HS-635HB', 6.0, 5.98, 49.9; ...
                    'HS-805MG', 6.0, 24.7, 197; ...
                    'HS-81', 6.0, 3.02, 16.4; ...
                    'HS-945MG', 6.0, 11.0, 55.9};
                    
torque_lowV = cell2mat(servo_data_lowV(:,3));
mass_lowV = cell2mat(servo_data_lowV(:,4));

torque_highV = cell2mat(servo_data_highV(:,3));
mass_highV = cell2mat(servo_data_highV(:,4));

% Case 1: V = 4.8V

coeffs_lowV = polyfit(torque_lowV, mass_lowV, 3);
torque_lowV_curve = linspace(min(torque_lowV), max(torque_lowV));
line_best_fit_lowV = polyval(coeffs_lowV, torque_lowV_curve);

% Sample Plot
% figure(1)
% scatter(torque_lowV, mass_lowV)
% hold on
% plot(torque_lowV_curve, line_best_fit_lowV)
% xlabel('Torque [kg-cm]'); ylabel('Mass [g]')
% title('Torque-Mass Trend: V=4.8V')

% Case 2: V = 6.0V

coeffs_highV = polyfit(torque_highV, mass_highV, 3);
torque_highV_curve = linspace(min(torque_highV), max(torque_highV));
line_best_fit_highV = polyval(coeffs_highV, torque_highV_curve);

% Sample Plot
% figure(2)
% scatter(torque_highV, mass_highV)
% hold on
% plot(torque_highV_curve, line_best_fit_highV)
% xlabel('Torque [kg-cm]'); ylabel('Mass [g]');
% title('Torque-Mass Trend: V=6.0V')

if servo_voltage == 4.8
    massServo = coeffs_lowV(4) + coeffs_lowV(3).*servo_torque + ...
        coeffs_lowV(2).*servo_torque.^2 + coeffs_lowV(1).*servo_torque.^3;
else
    massServo = coeffs_highV(4) + coeffs_highV(3).*servo_torque + ...
        coeffs_highV(2).*servo_torque.^2 + coeffs_highV(1).*servo_torque.^3;
end

massServo = massServo / 1000; % g --> kg
end
