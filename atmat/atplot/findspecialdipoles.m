function varargout = findspecialdipoles(dipoles,rok)
%FINDSPECIALDIPOLES identify bending magnet types
%[bendType, ColCode] = findspecialdipoles(dipoles,rok)
%
%  Helper function called by atplotsyn to identify different kind of
%  dipoles, i.e. AntiBends (AB), Transverse Gradient Bends (TGB) and
%  ordinary Bends (B). These are assigned different colour codes for an easy
%  visual check in atplot
% 
%  See also atplot, atplotsyn
%
% Author MA 13092023

indx = find(dipoles==1);
bendType = cell(numel(indx),1);
ColCode = cell(numel(indx),1);
for i = 1:length(indx)
    bendType{i} = 'B';
    ColCode{i}  = [0.5 0.5 1]; 
    if rok{indx(i)}.BendingAngle < 0 
        bendType{i} = 'AB';
        ColCode{i}  = [0.7 0.3 1]; 
    elseif rok{indx(i)}.PolynomB(2) ~= 0
        bendType{i} = 'TGB';
        ColCode{i}  = [0.8 0.8 1]; 
    end

end
if nargout > 0
    varargout{1} = bendType;
end

if nargout > 1
    varargout{2} = ColCode;
end
end