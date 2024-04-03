function Qx=nux(ring, periodicity)
% get value of horizontal dispersion for  Seq(indx)
% % global THERING
% % THERING = ring;

% % AAA=atsummary;
% % Qx=mod(AAA.tunes(1) * periodicity,1);
% % %Qy=AAA.tunes(2) * periodicity;

if nargin<2
    periodicity =1;
end

% if nux calculation is non physical, give it an absurd value
try
    [ringdata,lindata]=atlinopt6(ring,1:length(ring)+1); 
    Qx = lindata(end).mu(1)/2/pi * periodicity; % 
    if ~isreal(Qx)
        Qx = 1e19;
    end
catch
    Qx = 1e19;
end

