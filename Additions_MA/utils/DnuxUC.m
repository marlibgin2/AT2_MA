function DeltaNuxUC=DnuxUC(ring, periodicity)
% get value of horizontal dispersion for  Seq(indx)
% % global THERING
% % THERING = ring;

% % AAA=atsummary;
% % Qx=mod(AAA.tunes(1) * periodicity,1);
% % %Qy=AAA.tunes(2) * periodicity;

if nargin<2
    periodicity =1;
end
if 1 == 0 % scheme specific for medmax
S1hi = findcells(ring,'FamName','S1h'); 
%S1h_S = findspos(rrr, S1hi); 

[ringdata,lindata]=atlinopt6(ring,1:length(ring)+1);
delta = (lindata(S1hi(9)).mu - lindata(S1hi(8)).mu)/2/pi;
DeltaNuxUC = delta(1); 
end

S4i = findcells(ring,'FamName','S4'); 
S5i = findcells(ring,'FamName','S5'); 
[ringdata,lindata]=atlinopt6(ring,1:length(ring)+1);
delta = (lindata(S5i(2)).mu - lindata(S5i(1)).mu)/2/pi;
DeltaNuxUC = delta(1); 



end