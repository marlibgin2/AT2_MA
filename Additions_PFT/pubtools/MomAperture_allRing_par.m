function [map_l,map_h]=MomAperture_allRing_par(THERING,points,varargin)
% all Ring momentum aperture
% points=1:10:length(THERING);
map_h=zeros(length(points),1);
map_l=zeros(length(points),1);

%% Input argument parsing

deltalimit    = getoption(varargin,'deltalimit',0.1);
initcoord     = getoption(varargin,'initcoord',[1E-6 1E-6]);
delta         = getoption(varargin,'delta',0.0);
precdelta     = getoption(varargin,'precdelta',0.0);
deltastepsize = getoption(varargin,'deltastepsize',0.005);
splits        = getoption(varargin,'splits',2);
split_step_divisor = getoption(varargin,'split_step_divisor',10);
nturns        = getoption(varargin,'nturns',100);

%% Calculates Momentum Aperture
parfor i=1:length(points)
   % disp([i length(points) i/length(points)*100])
    %cycle ring
     THERING_cycl=[THERING(points(i):end); THERING(1:points(i)-1)]';
        try
            map_h(i)=momentum_aperture_at(THERING_cycl,deltalimit,initcoord,...
                delta,precdelta,deltastepsize,splits,split_step_divisor,nturns);
        catch
            map_h(i)=0;
        end
        
        try
            map_l(i)=momentum_aperture_at(THERING_cycl,-deltalimit,initcoord,...
                delta,precdelta,-deltastepsize,splits,split_step_divisor,nturns);
        catch
            map_l(i)=0;
        end
       
end

return