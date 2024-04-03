function [ex]=emix(ring)
% get value of horizontal dispersion for  Seq(indx)
global THERING
THERING = ring;

% if emittance is non-physical, give it a huge value ... 
try
    AAA=atsummaryMA;
    ex=AAA.naturalEmittance;
catch
    disp('exception in emix found ... ')
    ex=1e19;
end

