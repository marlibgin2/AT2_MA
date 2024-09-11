function [map_l,map_h]=calcLMA_raw(RING,points,varargin)
% Parallel computation of local momentum aperture over a lattice.
% This is intended as a lower level function that is called by the higher
% level wrapper function "calcLMA"
%
%% Inputs
% Mandatory input arguments
% RING : AT  lattice array
% Ipos: array of lattice positions at which LMA is to be calculated
% 
% Optional input parameters
% input with format ('paramater',value)
%
% deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
% initcoord: initial coordinates [x0 y0]
% delta: initial guess for momentum aperture 
% deltastepsize: step size for LMA search;
% splits : number of iterations of step division
% split_step_divisor: factor to reduce step size at each iteration
% nturns: number of turns. If nan then number of turns is chosen as 1.2/Qs
%                         this is handled by the momentum_aperture_at
%                         function
% verbose : defines level of verbose output, default=0, i.e. no output
%
%% Outputs
% map_l : negative LMA
% map_h : positive LMA
%% Usage examples
% [map_l,map_h] = calcLMA_raw(RING,Ipos,'nturns',500)
%

%% History
% PFT 2024/03/08. 
% PFT 2024/06/16: changed handling of verbose level
<<<<<<< HEAD

%% Input argument parsing

deltalimit    = getoption(varargin,'deltalimit',0.1);
initcoord     = getoption(varargin,'initcoord',[30E-6 30E-6]);
delta         = getoption(varargin,'delta',0.1);
deltastepsize = getoption(varargin,'deltastepsize',0.001);
splits        = getoption(varargin,'splits',2);
split_step_divisor = getoption(varargin,'split_step_divisor',10);
nturns        = getoption(varargin,'nturns',500);
verboselevel  = getoption(varargin,'verbose',0);
=======
% PFT 2024/07/30: changed handling of initcoord

%% Input argument parsing

deltalimit      = getoption(varargin,'deltalimit',0.1);
initcoord       = getoption(varargin,'initcoord',[]);
delta           = getoption(varargin,'delta',0.1);
deltastepsize   = getoption(varargin,'deltastepsize',0.001);
splits          = getoption(varargin,'splits',2);
split_step_divisor = getoption(varargin,'split_step_divisor',10);
nturns          = getoption(varargin,'nturns',1024);
verboselevel    = getoption(varargin,'verbose',0);
>>>>>>> MAXIV_addition

%% Calculates Momentum Aperture
map_h=zeros(length(points),1);
map_l=zeros(length(points),1);
<<<<<<< HEAD
=======

% Adjusts the initial phase to the synchronous phase if input is nan
if (not(isempty(initcoord)))
    if(isnan(initcoord(6)))
        initco=findorbit6(RING,1);
        initcoord(6)=initco(6);
    end
end

>>>>>>> MAXIV_addition
parfor i=1:length(points)
    %cycle ring
     RING_cycl=[RING(points(i):end); RING(1:points(i)-1)];
        try
            if (verboselevel>0)
                fprintf('%s calcLMA_raw: calculating positive LMA \n', datetime);
            end
            if (not(isnan(nturns)))
                map_h(i)=momentum_aperture_at(RING_cycl,deltalimit,initcoord,...
                    delta,0.0,deltastepsize,splits,split_step_divisor,nturns);
            else
                map_h(i)=momentum_aperture_at(RING_cycl,deltalimit,initcoord,...
                    delta,0.0,deltastepsize,splits,split_step_divisor);
            end
<<<<<<< HEAD
        catch
            map_h(i)=0;
=======
        catch ME
            map_h(i)=0;
            fprintf('%s calcLMA_raw: Error in momentum_aperture_at \n', datetime);
            fprintf('Error message was:%s \n',ME.message);
>>>>>>> MAXIV_addition
        end
        
        try
            if (verboselevel>0)
                 fprintf('%s calcLMA_raw: calculating negative LMA \n', datetime);
            end
            if (not(isnan(nturns)))
                map_l(i)=momentum_aperture_at(RING_cycl,-deltalimit,initcoord,...
                -delta,0.0,-deltastepsize,splits,split_step_divisor,nturns);
            else
                map_l(i)=momentum_aperture_at(RING_cycl,-deltalimit,initcoord,...
                -delta,0.0,-deltastepsize,splits,split_step_divisor);
            end
        catch
            map_l(i)=0;
<<<<<<< HEAD
=======
            fprintf('%s calcLMA_raw: Error in momentum_aperture_at \n', datetime);
            fprintf('Error message was:%s \n',ME.message);
>>>>>>> MAXIV_addition
        end
       
end

return