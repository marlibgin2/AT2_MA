function Qy=nuy(ring, periodicity)
% get value of horizontal dispersion for  Seq(indx)
% % % global THERING
% % % THERING = ring;

% % % AAA=atsummary;
% % % %Qx=AAA.tunes(1) * periodicity;
% % % Qy=mod(AAA.tunes(2) * periodicity,1);

if nargin<2
    periodicity =1;
end
% if nuy calculation is non physical, give it an absurd value
try
    [ringdata,lindata]=atlinopt6(ring,1:length(ring)+1);
    Qy = lindata(end).mu(2)/2/pi * periodicity; %
    if ~isreal(Qy)
        Qy = 1e19;
    end
catch
    Qy = 1e19;
end