function PlotTwiss(LATTICE,split)
%UNTITLED2 Plots beta and dispersion of periodic lattice 
%   splitting each element by factor "spli"
LAT=LATTICE;
if (split>1)
    for i=1:split
        LAT=atinsertelems(LAT,[1:length(LAT)],0.5,[]);
    end
end
lindata=atlinopt(LAT,0.0,1:length(LAT)+1);
PlotBetaDisp(lindata);
end

