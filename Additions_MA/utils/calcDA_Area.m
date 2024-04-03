function [Area, a] = calcDA_Area(DA)
%
% external fn 
% 

disp('... this is a calcDA_area external from the main MOGA')
O = [0, 0];
X = DA(:,1);
Y = DA(:,2);

A = [0, 0];
for i = 1:length(DA)-1
    B = DA(i,  :);
    C = DA(i+1,:);
    a(i) = (A(1) * (B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))/2;
end

Area = sum(a);