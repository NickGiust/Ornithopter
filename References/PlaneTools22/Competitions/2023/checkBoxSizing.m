function [plane, box] = checkBoxSizing(plane)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% INPUTS

maxDim = 62*0.0254;
box.thickness = 1/2 * 0.0254; % 0.5 in --> m


% antennaAdapterSpan = 0.035; % antenna adapter additional width (contributing to span)
antennaAdapterSpan = 0.0; % antenna adapter NOT ATTACHED TO WING while in box
% adding adapter width to span because adapter is preinstalled on all wingtips
outerSectionSpan = plane.outerSectionSpan + antennaAdapterSpan;
antennaAdditionalWidth = 1*0.0254;


%% PIZZA BOX WITH DISASSEMBLED (TWO-PIECE) TAIL

box.width = 2*plane.croot + plane.wFuse + 2*box.thickness;
box.height = max([plane.hFuse, 2*plane.croot*plane.tcRatio]) + 2*box.thickness; % square shape
availableInternalLength = maxDim - box.width - box.height - 2*box.thickness;

if max([outerSectionSpan, plane.lFuse]) <= availableInternalLength % wing sections will fit
    box.type = "PIZZA BOX w/ DISASSEMBLED TAIL";
    plane.planeFits_boolean = true;
    %disp("Pizza box plane fits")

    if availableInternalLength >= max([outerSectionSpan, plane.lFuse]) + antennaAdditionalWidth
        % give an inch of length to the width to store antenna diagonally
        box.height = box.height + antennaAdditionalWidth; % assume an extra inch of box width for antenna
        box.length = maxDim - box.width - box.height;
        available_length_antenna = sqrt((box.width-2*box.thickness)^2+(box.length-2*box.thickness)^2) - 0.02;
    else
        box.length = maxDim - box.width - box.height;
        available_length_antenna = box.length - 2*box.thickness;
        disp("Pizza box antenna does not fit")
    end
else % wing sections do not fit in box (IMPOSSIBLE - fail the plane)
    plane.planeFits_boolean = false;
    %disp("PLANE DOES NOT FIT IN SAFE BOX")
    box.length = max([outerSectionSpan plane.lengthPayload3]) + 2*box.thickness;
    available_length_antenna = 1; % doesn't matter b/c mission scores will be 0 if plane does not fit
end

box.totalDimensionMeters = box.length + box.width + box.height;
box.totalDimensionInches = box.totalDimensionMeters/0.0254;

if plane.lengthPayload3 <= available_length_antenna % antenna fits in box
    plane.antennaFits_boolean = true;
else
    plane.antennaFits_boolean = false;
    disp("ANTENNA DOES NOT FIT IN PIZZA BOX")
end






% 
% 
% %% PIZZA BOX WITH ASSEMBLED (SINGLE-PIECE) TAIL
% 
% if plane.planeFits_boolean == false || plane.antennaFits_boolean == false
%     % PIZZA BOX w/disassembled tail did not work, try PIZZA BOX w/assembled tail
%     
%     box.width = 2*plane.croot + plane.wFuse + 2*box.thickness;
%     box.height = max([plane.hFuse, 2*plane.croot*plane.tcRatio, plane.ch]) + 2*box.thickness; % square shape
%     availableInternalLength = maxDim - box.width - box.height - 2*box.thickness;
%     
%     if outerSectionSpan <= availableInternalLength % wing sections will fit
%         box.type = "PIZZA BOX w/ ASSEMBLED TAIL";
%         plane.planeFits_boolean = true;
%         disp("Pizza box plane fits")
%     
%         if availableInternalLength >= outerSectionSpan + antennaAdditionalWidth
%             % give an inch of length to the width to store antenna diagonally
%             box.height = box.height + antennaAdditionalWidth; % assume an extra inch of box width for antenna
%             box.length = maxDim - box.width - box.height;
%             available_length_antenna = sqrt((box.width-2*box.thickness)^2+(box.length-2*box.thickness)^2);
%         else
%             box.length = maxDim - box.width - box.height;
%             available_length_antenna = box.length - 2*box.thickness;
%             disp("Pizza box antenna does not fit")
%         end
%     else % wing sections do not fit in box (IMPOSSIBLE - fail the plane)
%         plane.planeFits_boolean = false;
%         %disp("PLANE DOES NOT FIT IN SAFE BOX")
%         box.length = max([outerSectionSpan plane.lengthPayload3]) + 2*box.thickness;
%         available_length_antenna = 1; % doesn't matter b/c mission scores will be 0 if plane does not fit
%     end
%     
%     box.totalDimensionMeters = box.length + box.width + box.height;
%     box.totalDimensionInches = box.totalDimensionMeters/0.0254;
%     
%     if plane.lengthPayload3 <= available_length_antenna % antenna fits in box
%         plane.antennaFits_boolean = true;
%     else
%         plane.antennaFits_boolean = false;
%         disp("ANTENNA DOES NOT FIT IN PIZZA BOX")
%     end
% 
% end

%{
%% WHIRLPOOL BOX DESIGN

if plane.planeFits_boolean == false || plane.antennaFits_boolean == false % if PIZZA BOX did not work, try WHIRLPOOL design

    box.width = plane.croot + (1.5*plane.tcRatio*plane.croot) + 2*box.thickness;
    box.height = box.width; % square shape
    availableInternalLength = maxDim - box.width - box.height - 2*box.thickness;
    
    if outerSectionSpan <= availableInternalLength % wing sections will fit
        box.type = "WHIRLPOOL";
        plane.planeFits_boolean = true;
        % availableInternalLength = max([plane.b*outboardSecFrac plane.lengthPayload3]);
        % box.length = availableInternalLength + 2*box.thickness;
        
        box.length = maxDim - box.width - box.height;
        available_length_antenna = box.length - 2*box.thickness;
    else % wing sections do not fit in box (IMPOSSIBLE - fail the plane)
        plane.planeFits_boolean = false;
        % disp("PLANE DOES NOT FIT IN WHIRLPOOL")
        box.length = max([outerSectionSpan plane.lengthPayload3]) + 2*box.thickness;
        available_length_antenna = 1; % doesn't matter b/c mission scores will be 0 if plane does not fit
    end
    
    box.totalDimensionMeters = box.length + box.width + box.height;
    box.totalDimensionInches = box.totalDimensionMeters/0.0254;
    
    if plane.lengthPayload3 <= available_length_antenna % antenna fits in box
        plane.antennaFits_boolean = true;
    else
        plane.antennaFits_boolean = false;
        %disp("ANTENNA DOES NOT FIT IN WHIRLPOOL")
    end

end

%% SAFE BOX DESIGN

if plane.planeFits_boolean == false || plane.antennaFits_boolean == false % if WHIRLPOOL did not work, try SAFE BOX design

    box.width = max([plane.croot, max([plane.wFuse plane.hFuse])]) + 2*box.thickness;
    % fuselage stored length-wise in the box, center wing section span
    % contained within box width
    box.height = (4*plane.tcRatio*plane.croot) + min([plane.wFuse plane.hFuse]) + 2*box.thickness;
    % 4x wing thicknesses due to four outer wing sections
    availableInternalLength = maxDim - box.width - box.height - 2*box.thickness;
    
    if outerSectionSpan <= availableInternalLength % wing sections will fit
        box.type = "SAFE";
        plane.planeFits_boolean = true;
        % availableInternalLength = max([plane.b*outboardSecFrac plane.lengthPayload3]);
        % box.length = availableInternalLength + 2*box.thickness;
        if availableInternalLength >= outerSectionSpan + antennaAdditionalWidth
            % give an inch of length to the width to store antenna diagonally
            box.width = box.width + antennaAdditionalWidth; % assume an extra inch of box width for antenna
            box.length = maxDim - box.width - box.height;
            available_length_antenna = sqrt((box.height-2*box.thickness)^2+(box.length-2*box.thickness)^2);
        else
            box.length = maxDim - box.width - box.height;
            available_length_antenna = box.length - 2*box.thickness;
        end
    else % wing sections do not fit in box (IMPOSSIBLE - fail the plane)
        plane.planeFits_boolean = false;
        %disp("PLANE DOES NOT FIT IN SAFE BOX")
        box.length = max([outerSectionSpan plane.lengthPayload3]) + 2*box.thickness;
        available_length_antenna = 1; % doesn't matter b/c mission scores will be 0 if plane does not fit
    end
    
    box.totalDimensionMeters = box.length + box.width + box.height;
    box.totalDimensionInches = box.totalDimensionMeters/0.0254;
    
    if plane.lengthPayload3 <= available_length_antenna % antenna fits in box
        plane.antennaFits_boolean = true;
        if plane.planeFits_boolean == true
            disp("PLANE FIT IN SAFE BOX")
        end

    else
        plane.antennaFits_boolean = false;
        %fprintf("ANTENNA DOES NOT FIT IN SAFE BOX\n")
    end

else

    disp("PLANE FIT IN WHIRLPOOL")

end

%}