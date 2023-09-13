function [bendType, ColCode] = FindSpecialDipoles(dipoles,rok)
% -----------------------------------------------------------------------
% function called by atplotsyn to identify differnet kind of dipoles (AB,
% TGB, B) and assign different color codes for an easy visual check in
% atplot
% MA: 13092023
% -----------------------------------------------------------------------
indx = find(dipoles==1)
for i = 1:length(indx)
    bendType{i} = 'B';
    ColCode{i}  = [0.5 0.5 1]; 
    if rok{indx(i)}.BendingAngle < 0 ;
        bendType{i} = 'AB';
        ColCode{i}  = [0.7 0.3 1]; 
    elseif rok{indx(i)}.PolynomB(2) ~= 0
        bendType{i} = 'TGB';
        ColCode{i}  = [0.8 0.8 1]; 
    end

end

end