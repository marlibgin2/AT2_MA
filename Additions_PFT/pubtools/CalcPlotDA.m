function DA=CalcPlotDA(varargin)
%
% Calculates and plots (if plotf=1) DA for Lattice "RING"
% Parameters for DA calculation are all specified in structure
% LatticeOptDat.DAoptions
% Usage: DA=CalcPlotDA(RING,DAoptions,plotf)
%
%
%% Input argument parsing
[RING,DAoptions]=getargs(varargin,[],[]);
plotf          = any(strcmpi(varargin,'plot'));
dp             = getoption(varargin,'dp',DAoptions.dp);
z0             = getoption(varargin,'z0',DAoptions.z0);

DAoptions.dp=dp;
DAoptions.z0=z0;
% Parameters for dynamic aperture calculation
%

%% Checks backward compatibilty
if (~isfield(DAoptions,'xmaxdas'))
    DAoptions.xmaxdas = 0.007;
end

if (~isfield(DAoptions,'xmindas'))
    DAoptions.xmindas = -0.015;
end

if (~isfield(DAoptions,'ymaxdas'))
    DAoptions.ymaxdas = 0.003;
end
   
if(~(isfield(DAoptions,'XmaxDA')))
    DAoptions.XmaxDA = 0.015;
end

if(~(isfield(DAoptions,'YmaxDA')))
    DAoptions.YmaxDA     = 0.004;
end

XmaxDA = DAoptions.XmaxDA;
YmaxDA = DAoptions.YmaxDA;

if(~(isfield(DAoptions,'DAmode')))
    DAoptions.DAmode = 'border';
end
DAmode = DAoptions.DAmode;

if (strcmp(DAmode,'grid'))
%
%% Recalculates X0da and Y0da in case the data in DAoptions 
% is not consistent (e.g. is not the same as ina previous MOGA run)
%
    npdax      = DAoptions.npdax;   % number of grid points in x direction is 2*npdax+1
    npday      = DAoptions.npday;   % number of grid points in y direction is  npday+1
    npDA       = DAoptions.npDA;    % total numbr of grid points
    dx = XmaxDA/npdax; % grid stepsize in x [m]
    dy = YmaxDA/npday;   % grid stepsize in y [m]
    dxdy = dx*dy;
    X0da = zeros(npDA,1);  % horizontal coordinates of grid points [m]
    Y0da = zeros(npDA,1);  % vertical coordinates of grid points [m]
    k= 1;
    for i=0:npday 
        for j=0:2*npdax
        X0da(k) = -XmaxDA+dx*j;
        Y0da(k) =  dy*i;
        k=k+1;
        end
    end
    DAoptions.X0da=X0da;
    DAoptions.Y0da=Y0da;
    DAoptions.dxdy=dxdy;
end
PC=load('PC.mat');      %to prevent matlab from complaining about variable name being the same as script name.
PhysConst = PC.PC;      %Load physical constants

%% Calculates and plots DA
angle=atgetfieldvalues(RING,'BendingAngle');
isdipole=isfinite(angle) & (angle~=0);

try
   rpara = atsummary_fast(RING,isdipole);
   etax = rpara.etax;
   if (isnan(z0))
       ats=atsummary(RING);
       z0 = PhysConst.c*(ats.syncphase-pi)/(2*pi*ats.revFreq*ats.harmon); %if z0 not givem choose the synchronous phase
       DAoptions.z0=z0;
   end
   [DA,DAV] = calcDA_fast(RING,DAoptions,etax,rpara.beta0(1),rpara.beta0(2));
   if (plotf)
     switch DAmode
       case 'border'
           figure;plot(DAV(:,1)*1000,DAV(:,2)*1000,'-ob');
           xlabel('X [mm]'); ylabel('Y [mm]');grid;
           xlim([-XmaxDA XmaxDA]*1000);ylim([0 YmaxDA]*1000);

         case 'grid'
           DAM = zeros(npday+1,2*npdax+1);
           k= 1;
           for i=0:npday
              for j= 1:2*npdax+1
                DAM(npday+1-i,j)=DAV(k);
                k=k+1;
               end
            end
            DAM=DAM*255;
            map=[0 0.75 0; 1 1 1];
            figure;image([-XmaxDA*1000,XmaxDA*1000],[YmaxDA*1000,0],DAM);
            ax=gca;
            ax.YDir='normal';
            colormap(map);
            xlabel('X[mm]');
            ylabel('Y[mm]');    
            xlim([-XmaxDA XmaxDA]*1000);ylim([0 YmaxDA]*1000);grid;  
     end
   end
catch ME
     fprintf('Error calculating Dynamic Aperture \n');
     fprintf('Error message was:%s \n',ME.message);
     DA=NaN;
end