function [massMotorMount] = calcMass_MotorMount(plane)
%Lets weigh the motor mount baby!!

% Constants
t = 0.25 * 0.0254; % [in -> m]
density_ply = 680; %kg/m3
epoxyFrac = 1.3; % additional mass from epoxy
numMotors = plane.nMotors;
x = plane.motorLength*1.1; % [m]
d = plane.motorDiam*1.1; % [m]

if rem(numMotors,2) == 1
    numMotors_Fuse = 1;
    numMotors_Wings = numMotors - 1;
else
    numMotors_Fuse = 0;
    numMotors_Wings = numMotors;
end

%% Wing Motor Mount Mass
if numMotors_Wings ~= 0
%     SA = 4*(x+t)*(d+t) + (d+t)^2; % [m^2]
    SA = d^2; % [m^2] Back Mounted
    massMotorMount_Wing = numMotors_Wings*(SA*t*density_ply)*(epoxyFrac)*(0.75); % [kg]
else
    massMotorMount_Wing = 0;
end

%% Fuselage Motor Mount Mass
build_method = plane.fuseBuildMethod;

if numMotors_Fuse ~= 0
    if strcmpi(build_method, 'Builtup')
        Hl = input('Height of Longeron: '); % [m]
        Hf = plane.hFuse; % [m]
        Nl = 1.3*Hf; % [m]
        massMotorMount_Fuse = (2*Hl*Nl*t*density_ply)*epoxyFrac; % [kg]
    elseif strcmpi(build_method, 'Monocoque')
        massMotorMount_Fuse = (d^2)*(0.125*0.0254)*density_ply; % [kg]
    end
else
    massMotorMount_Fuse = 0;
end
 
%% Mass Total

massMotorMount = numMotors_Fuse*massMotorMount_Fuse + numMotors_Wings*massMotorMount_Wing;

end