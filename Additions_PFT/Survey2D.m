function [x2d, y2d, a2d] = Survey2D(LATTICE,STARTANGLE)
% Determine 2-d geometry of the LATTICE
%
% Mod by Pedro F. Tavares
% 2023/02/05
% to allow for an open lattice, that does not closes on itself, e.g. a 
% single cell or single achromat
%
NumElements = length(LATTICE);
x2d = zeros(1,NumElements+1);
y2d = zeros(1,NumElements+1);
a2d = zeros(1,NumElements+1); % angle of orbit in radians
a2d(1) = STARTANGLE;
%for en = 1:NumElements-1
 for en = 1:NumElements
    if isfield(LATTICE{en},'BendingAngle')
        ba = LATTICE{en}.BendingAngle; % bending angle in radians
    else
        ba = 0;
    end

    if ba == 0
        Lt = LATTICE{en}.Length;
        Lp = 0;
    else
        Lt = LATTICE{en}.Length*sin(ba)/ba;
        Lp = -LATTICE{en}.Length*(1-cos(ba))/ba;
    end

    x2d(en+1) = x2d(en) + Lt*cos(a2d(en)) - Lp*sin(a2d(en));
    y2d(en+1) = y2d(en) + Lt*sin(a2d(en)) + Lp*cos(a2d(en));
    a2d(en+1)=a2d(en) - ba;

end
