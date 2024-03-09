function [map_l,map_h]=calcLMA_raw(RING,points,varargin)
% Parallel computation of local momentum aperture over a lattice.
% This is intended as a lower level function that is called by the higher
% level wrapper function "calcLMA"
%
%% Usage example
% [map_l,map_h] = calcLMA_raw(RING,Ipos,'nturns',500)
%
%% Mandatory input arguments
% RING : AT  lattice array
% Ipos: array of lattice positions at which LMA is to be calculated
% 
%% Optional input parameters
% input with format ('paramater',value)
%
% deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
% initcoord: initial coordinates [x0 x0p y0 x0p delta z0]'
% delta: initial guess for momentum aperture 
% deltastepsize: step size for LMA sereach;
% splits : number of iterations of step division
% split_step_divisor: factor to reduce step size at each iteration
% nturns: number of turns. If nan then number of turns is chosen as 1.2/Qs
%                         this is handled by the momentum_aperture_at
%                         function
%% Output parameters
% map_l : negative LMA
% map_h : positive LMA

map_h=zeros(length(points),1);
map_l=zeros(length(points),1);

%% Input argument parsing

deltalimit    = getoption(varargin,'deltalimit',0.1);
initcoord     = getoption(varargin,'initcoord',[1E-6 1E-6]);
delta         = getoption(varargin,'delta',0.1);
deltastepsize = getoption(varargin,'deltastepsize',0.001);
splits        = getoption(varargin,'splits',2);
split_step_divisor = getoption(varargin,'split_step_divisor',10);
nturns        = getoption(varargin,'nturns',500);

%% Calculates Momentum Aperture
parfor i=1:length(points)
    %cycle ring
     RING_cycl=[RING(points(i):end); RING(1:points(i)-1)];
        try
            if (not(isnan(nturns)))
                map_h(i)=momentum_aperture_at(RING_cycl,deltalimit,initcoord,...
                    delta,0.0,deltastepsize,splits,split_step_divisor,nturns);
            else
                map_h(i)=momentum_aperture_at(RING_cycl,deltalimit,initcoord,...
                    delta,0.0,deltastepsize,splits,split_step_divisor);
            end
        catch
            map_h(i)=0;
        end
        
        try
            if (not(isnan(nturns)))
                map_l(i)=momentum_aperture_at(RING_cycl,-deltalimit,initcoord,...
                -delta,0.0,-deltastepsize,splits,split_step_divisor,nturns);
            else
                map_l(i)=momentum_aperture_at(RING_cycl,-deltalimit,initcoord,...
                -delta,0.0,-deltastepsize,splits,split_step_divisor);
            end
        catch
            map_l(i)=0;
        end
       
end

return