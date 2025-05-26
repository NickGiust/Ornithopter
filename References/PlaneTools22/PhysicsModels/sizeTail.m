function [plane] = sizeTail(plane)
    % Arbitrary moment arm
    plane.lh = 0.6*plane.b; % old value
%     plane.lh = 0.6*plane.lFuse; % @DAVID, does this make sense?

    % Horizontal sizing
    plane.Sh = plane.Vh*plane.S*plane.c/plane.lh;
    plane.bh = sqrt(plane.ARh*plane.Sh);
    plane.ch = plane.Sh./plane.bh;

    % Vertical sizing
    plane.lv = plane.lh;
    plane.Sv = plane.Vv*plane.S*plane.b/plane.lv;
    plane.cv = plane.ch;
    plane.bv = plane.Sv/plane.cv;
    plane.ARv = plane.bv.^2./plane.Sv;

    % Tail parameters (chord and span from hor, height from vert)
    plane.cTail = plane.ch;
    plane.bTail = plane.bh;
    plane.hTail = plane.bv;
end