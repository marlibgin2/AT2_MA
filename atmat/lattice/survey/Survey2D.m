function [s2d,x2d, y2d, a2d,baa, ban] = Survey2D(LATTICE,STARTANGLE)
% Determine 2-d geometry of the LATTICE
%
% Mod. by Pedro F. Tavares
% 2023/02/05
% to allow for an open lattice, that does not close on itself, e.g. a 
% single cell or single achromat
%
% outputs:
% s2d: orbit length
% x2d: horizontal coordinate [m]
% y2d: vertical coordinate [m]
% a2d: angle [rad]
% baa: change in angle [rad]
% ban: Element family if bending angle is not zero
%
NumElements = length(LATTICE);
x2d = zeros(1,NumElements+1);
y2d = zeros(1,NumElements+1);
a2d = zeros(1,NumElements+1); % angle of orbit in radians
baa = zeros(1,NumElements);
s2d = findspos(LATTICE,1:NumElements+1);
ban = cell(NumElements);
a2d(1) = STARTANGLE;
for en = 1:NumElements
    if isfield(LATTICE{en},'BendingAngle')
        ba = LATTICE{en}.BendingAngle; % bending angle in radians
        ban{en}=LATTICE{en}.FamName;
    else
        ba = 0;
        ban{en}='NF';
    end

    if ba == 0
        Lt = LATTICE{en}.Length;
        Lp = 0;
    else
        Lt = LATTICE{en}.Length*sin(ba)/ba;
        Lp = -LATTICE{en}.Length*(1-cos(ba))/ba;
    end
    baa(en)=ba;
    x2d(en+1) = x2d(en) + Lt*cos(a2d(en)) - Lp*sin(a2d(en));
    y2d(en+1) = y2d(en) + Lt*sin(a2d(en)) + Lp*cos(a2d(en));
    a2d(en+1) = a2d(en) - ba;
end
