function [dx,dxp,dy,dyp]=getDispersion2(RING,bpmindx)
% function [dx,dy]=getDispersion(THERING,bpmindx)
% 
% determines dispersion using twissring
% 

%[TD, tune, chrom] = twissring(RING,0,1:(length(RING)+1),'chrom',0.0001);
[ringdata,lindata]=atlinopt6(RING,1:length(RING)+1);
BETA = cat(1,lindata.beta);
S    = cat(1,lindata.SPos);
ETA  = cat(2,lindata.Dispersion);

dx  = ETA(1,bpmindx);
dxp = ETA(2,bpmindx);
dy  = ETA(3,bpmindx);
dyp = ETA(4,bpmindx);

return