function DAdist = calcDAdist(varargin)
% Calculates and plots Dynamic Aperture of a number of 
% lattice variants that differ only through the application of a given 
% error model (and possibly corresponding corrections).
% Tracking can be 6d or 4d as defined by the input lattice. 
% This is a higher level wrapper function
% that in turn calls the lower level function "calcDA_raw"
% 
%% Inputs
% Mandatory arguments
% RING : AT2 lattice array
% ErrorModel: structure generated by the errormodel function : 
%              if empty default is 
%              errormodel_DDRchallenging('gdran',0.0,'magalran',1.0,...
%                                'mulsys',1.0, 'mulran',1.0,...
%                                'bpmran',1.0, 'strran',1.0);
%
% DAoptions: Structure containing the following fields:
%     DAoptions.DAmode   = 'grid','border', 'smart_in' or 'smart_out'(plots are only produced if
%                           DAmode is NOT 'grid')
%      DAoptions.nturns   : number of turns, if nan use nsynchT synchrotron periods
%      DAoptions.nsyncT   : number of synchrotron periods to be used if nturns is n
%      DAoptions.betax0   : horizontal beta for normalization - if NaN, no normalization is done
%      DAoptions.betay0   : vertical beta for normalization - if NaN no normalization is done
%      DAoptions.xmaxdas  : limits of the range in which the DA border is searched
%      DAoptions.xmindas  : limits of the range in which the DA border is searched
%      DAoptions.ymaxdas  : limits of the range in which the DA border is searched
%      DAoptions.npd      : number of points along the energy axis for xdp and ydp 
%                           calculation modes
%      DAoptions.dpmin    : minimum energy deviation
%      DAoptions.dpmax    : maximum energy deviation
%
% Parameters for "border" DA calculation mode
%      DAoptions.r0      : initial guess [m]
%      DAoptions.nang    : number of angular steps
%      DAoptions.dp      : initial dp/p (6d tracking) or fixed dp/p (4d tracking)
%      DAoptions.z0      : initial longitudinal coordinate (6d tracking). Nan uses synchrnous phase
%      DAoptions.res     : resolution [m]
%      DAoptions.alpha   : da enlargement factor for DA border search
%
% Parameters for "grid" DA calculation mode
%      DAoptions.XmaxDA  : Horizontal range is -Xmax to Xmax [m]
%      DAoptions.YmaxDA  : Vertical range is 0 to Ymax [m]
%      DAoptions.npdax   : number of grid points in x direction is 2*npdax+1
%      DAoptions.npday   : number of grid points in y direction is  npday+
%
%            If DAoptions = [], hard-coded defaults are used.
%            Values of DAoptions fields are overridden if given explicitly 
%            as input in the form ('parameter', value)
%
% Optional arguments
%   all fields in DAoptions
%
%   ERlat: structure generated by the "generate_errlatt" function. If not
%   available or empty this structure is generated in this function by a
%   call to "generate_errlatt".
%
%   desc: descriptive string
%   mode : calculation mode, default = 'xy'. Possible values are
%       'xy'  : calculates DA in the xy plane for a given initial energy
%               deviation. Tracking may be 4d + energy or 6d as defined by
%               the input lattice.
%       'xydp':  calculates DA in the xyp planes. Tracking may be 4d + 
%               energy or 6d as defined by the input lattice.
%
%   corrorb: if true, perform orbit correction
%   corrtun: if true, perform tune correction
%   useORM0   : if true, sets the orbit correction to use the orbit respose
%               matrix for the unperturbed ring for all iterations, 
%               default=true
%
%   nseeds  : number of seeds, default = 10
%   tunfams : list of magnet families used for ring tun matching, default = {'Q1_b3','Q2_b3'}
%   nittune : number of iterations for tune matching, default = 10
%   TolTune : tolerance for tune matching, default = 1E-3
%   frac    : fraction for quad change in each tune fit iteration, defaut = 1.0
%
%   verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags
% plot : plots DA distribution;
% fulloutput : includes detailed data all data on each one of seeds.
%
%% Outputs
% DAdist structure with fields
%   DAdist.inputs echoes the inputs
%   DAdist.inputs.RING : input ring array
%   DAdist.inputs.ErrorModel 
%   DAdist.inputs.mode
%   DAdist.inputs.nseeds 
%   DAdist.inputs.tunfams 
%   DAdist.inputs.nittune 
%   DAdist.inputs.TolTune 
%   DAdist.inputs.frac
%   DAdist.inputs.corrorbf : correct orbit flag
%   DAdist.inputs.corrtunf : correct tune flag
%
%   DAdist.outputs.desc : : datetime + input description
%   DAdist.outputs.DAoptions : DAoptions actually used for the calculations
%   DAdist.outputs.DAs: (1Xnseeds+1) array of dynamic apertures for all seeds. 
%                     First point is the unperturbed lattice [mm**2]
%   DAdist.outputs.DAav: Average Dynamic aperture [mm**2]
%   DAdist.outputs.DAstd: DAStandard Deviation of Dynamics Apertures [mm**2]
%   DAdist.outputs.DAVs=DAVs: (nang+1X2*(nseeds+1)) array of dynamic aperture 
%                           border coordinates for all seeds
%                           (only for DAMode='border')
%   DAdist.outputs.orb0_stds: (6Xnseeds+1) array of closed orbit standard
%                           deviations for perturbed lattices before 
%                           correction. First point is the unperturbed
%                           lattice.
%   DAdist.outputsorb_stds  = (6Xnseeds+1) array of close orbit standard
%                           deviations for perturbed lattices after 
%                           correction. First point is the unperturbed
%                           lattice.
% Note: the outptus (RINGe,raparae,Itunese and Ftunese) below are empty 
%       unless the 'fulloutput' option is on
%   DAdist.outputs.RINGe:  (nseeds+1Xsize of RING) cell array of perturbed 
%                        lattices after correction. first is the
%                        unperturbed lattice.
%   DAdist.outputs.rparae: (nseeds+1X1) cell array of atsummaries for
%                        perturbed lattices after correction. 
%                        First is the unperturbed lattice.
%   DAdist.outputs.Itunese:(nseeds+1X1) cell array of tunes for the 
%                         perturbed lattices before correction
%                         First is the unperturbed lattice.
%   DAdist.outputs.Ftunese: (nseeds+1X1) cell array of tunes for the 
%                         perturbed lattices after correction
%                         First is the unperturbed lattice.
%   DAdist.outputs.stab   : (1xnseeds+1) array of integers. If one, seed is
%                           stable, if zero it is unstable
%   DAdist.outputs.survivalrate : pecentage of perturbed lattices that are
%                                 stable
%
%   DAdist.outputs.dps    : (npd X 1) array of momentum deviations
%   DAdist.outputs.DAxdppav:(npd X 1) array of average positive horizontal
%                          DA [m]
%   DAdist.outputs.DAxdppst:(npd X 1) array of standard deviation of 
%                          positive horizontal DA [m]
%   DAdist.outputs.DAxdpmav:(npd x 1) array of average negative horizontal
%                          DA [m]
%   DAdist.outputs.DAxdpmst:(npd x 1) array of standard deviation of 
%                          negative horizontal DA [m]
%   DAdist.outputs.DAydpav: (npd x 1) array of average vertical DA [m]
%   DAdist.outputs.DAydpst: (npd x 1) array of standard devitation of 
%                         average vertical DA [m]
%   DAdist.outputs.DAxdpsp = (npd X nseeds+1) array of positive horizontal
%                           DA [m]
%   DAdist.outputs.DAxdpsm = (npd X nseeds+1) array of negative horizontal
%                           DA [m]
%   DAdist.outputs.DAydps  = (npd X nseeds+1) array of vertical DA [m]
%
%   DAdist.outputs.telapsed: calculation time [s]
%
%% Usage examples
% DAdist = calcDAdist(RING,ErrorModel,DAoptions,'plot','verbose', 1);
% DAdist = calcDAdist(RING,ErrorModel,[],'nturns',1024);
% calcDAdist(RING,ErrorModel,[],'nturns',1024,'nseeds',10,'plot','corrorb',false);
% DAdist = calcDAdist(RING,ErrorModel,[],'nturns',810,'tunfams',{'Q1','Q2'},'frac',0.5);
% DAdist = calcDAdist(RING,ErrorModel,[],'mode,'xydp','nturns',810,'tunfams',{'Q1','Q2'});
% DAdist = calcDAdist(RING,ErrorModel,[],'mode,'xydp','nturns',810,'tunfams',{'Q1','Q2'},'Erlat',Erlat);

%% History
% PFT 2024/03/12
% PFT 2024/03/29: changed oputput to a structure, echoing also input
%                  parameters
% PFT 2024/04/23: added capabilty to fit tunes
% PFT 2024/04/27: added recording&plotting of rms orbit before and after correction
% PFT 2024/05/24: added updated fittuneRS parameters, DA calculation for
%                 different energies and more detailed info
%                 on perturbed lattices in the output structure
% PFT 2024/05/26: added more detailed info for xydp mode in output structure
% PFT 2024/06/04: added "smart" DA calculation mode
% PFT 2024/06/08: fixed description of verboselevel input
% PFT 2024/07/23: improved handling of unstable lattices (removed from mean
%                  values)
% PFT 2024/07/24: removed ring without errors from calculation of average
%                 and standard deviation
% PFT 2024/07/25: replaced code for error generation and correction by
%                 call to "generate_errlatt".
%                 added possibility of input of ERlat structure
%                 added possibility of fixing the response matrix for all
%                 seeds.
% PFT 2024/07/28: adapted to run mode "smart_in"
% PFT 2024/07/30: added handling of nan as input value for DAoptions.nturns
% PFT 2024/08/05: fixed bug - incorrect initilization of output vectors
%                 if the number of seeds was larger than the default (10)
% PFT 2024/08/07: fixed bug handling of nturns=nan
% PFT 2024/09/06: added handling of DAoptions.nsyncT

%% Input argument parsing
[RING,ErrorModel,DAoptions] = getargs(varargin,[],[],[]);
if (isempty(ErrorModel))
    ErrorModel=errormodel_DDRchallenging('gdran',1.0,'magalran',1.0,...
                                         'mulsys',1.0, 'mulran',1.0,...
                                         'bpmran',1.0, 'strran',1.0);
end
ERlat =  getoption(varargin,'ERlat',struct());

if (isempty(DAoptions))
    DAoptions.dp=0.0;
    DAoptions.z0=nan;
    DAoptions.DAmode='border';
    DAoptions.nturns=nan;
    DAoptions.nsyncT=3;
    DAoptions.betax0=nan;
    DAoptions.betay0=nan;
    DAoptions.xmaxdas=0.007;
    DAoptions.xmindas =-0.0150;
    DAoptions.ymaxdas = 0.006;
    DAoptions.npd=11;
    DAoptions.dpmin=-0.04;
    DAoptions.dpmax= 0.04;
    DAoptions.XmaxDA = 0.015;
    DAoptions.YmaxDA = 0.007;
    DAoptions.npdax = 64;
    DAoptions.npday = 64;
    DAoptions.r0 = 0.020;
    DAoptions.nang = 20;
    DAoptions.res = 5E-4; 
    DAoptions.alpha = 1.1;
    DAoptions.dx=nan;
    DAoptions.dy=nan;
    DAoptions.dxdy=nan;
    DAoptions.npDA = (2*DAoptions.npdax+1)*(DAoptions.npday+1); %total number of grid points
    DAoptions.X0da = zeros(DAoptions.npDA,1);  % horizontal coordinates of grid points [m]
    DAoptions.Y0da = zeros(DAoptions.npDA,1);  % vertical coordinates of grid points [m]
end

plotf            = any(strcmpi(varargin,'plot'));
plotorbrmsf      = any(strcmpi(varargin,'plotorbrms'));
fulloutputf      = any(strcmpi(varargin,'fulloutput'));
corrorbf         = getoption(varargin,'corrorb',true);
corrtunf         = getoption(varargin,'corrtun',true);
desc             = getoption(varargin,'desc','calcDAdist:');
useORM0f         = getoption(varargin,'useORM0',true);
verboselevel     = getoption(varargin,'verbose',0);
mode             = getoption(varargin,'mode','xy');
dp               = getoption(varargin,'dp',DAoptions.dp);
z0               = getoption(varargin,'z0',DAoptions.z0);
DAmode           = getoption(varargin,'DAmode',DAoptions.DAmode);
nturns           = getoption(varargin,'nturns',DAoptions.nturns);
if isfield(DAoptions,'nsyncT')
    nsyncT       = getoption(varargin,'nsyncT',DAoptions.nsyncT);
else
    nsyncT       = 3;
end
betax0           = getoption(varargin,'betax0',DAoptions.betax0);
betay0           = getoption(varargin,'betay0',DAoptions.betay0);
xmaxdas          = getoption(varargin,'xmaxdas',DAoptions.xmaxdas);
xmindas          = getoption(varargin,'xmindas',DAoptions.xmindas);
ymaxdas          = getoption(varargin,'ymaxdas',DAoptions.ymaxdas);
npd              = getoption(varargin,'npd',DAoptions.npd);
dpmin            = getoption(varargin,'dpmin',DAoptions.dpmin);
dpmax            = getoption(varargin,'dpmax',DAoptions.dpmax);
XmaxDA           = getoption(varargin,'XmaxDA',DAoptions.XmaxDA);
YmaxDA           = getoption(varargin,'YmaxDA',DAoptions.YmaxDA);
npdax            = getoption(varargin,'npdax',DAoptions.npdax);
npday            = getoption(varargin,'npday',DAoptions.npday);
r0               = getoption(varargin,'r0', DAoptions.r0);
nang             = getoption(varargin,'nang',DAoptions.nang);
res              = getoption(varargin,'res',DAoptions.res); 
alpha            = getoption(varargin,'alpha',DAoptions.alpha);

nseeds           = getoption(varargin,'nseeds',10);
tunfams          = getoption(varargin,'tunfams',{'Q1_b3','Q2_b3'});
nittune          = getoption(varargin,'nittune',10); 
TolTune          = getoption(varargin,'TolTune',1E-3); 
frac             = getoption(varargin,'frac',1.0); 


DAoptions.dp=dp;
DAoptions.z0=z0;
if (strcmpi(mode,'xydp'))
    DAmode='border';
    nang=2;
end
DAoptions.DAmode=DAmode;
DAoptions.nturns=nturns;
DAoptions.nsyncT=nsyncT;
DAoptions.betax0=betax0;
DAoptions.betay0=betay0;
DAoptions.xmaxdas=xmaxdas;
DAoptions.xmindas=xmindas;
DAoptions.ymaxdas=ymaxdas;
DAoptions.npd=npd;
DAoptions.dpmin=dpmin;
DAoptions.dpmax=dpmax;
DAoptions.XmaxDA=XmaxDA;
DAoptions.YmaxDA=YmaxDA;
DAoptions.npdax=npdax;
DAoptions.npday=npday;
DAoptions.r0=r0;
DAoptions.nang=nang;
DAoptions.res=res;
DAoptions.alpha=alpha;

dps=linspace(dpmin,dpmax,npd);

%
%% Recalculates X0da and Y0da in case the data in DAoptions 
% is not consistent (e.g. is not the same as in a previous MOGA run)
%
if (strcmp(DAmode,'grid'))
    npdax      = DAoptions.npdax;   % number of grid points in x direction is 2*npdax+1
    npday      = DAoptions.npday;   % number of grid points in y direction is  npday+1
    npDA       = DAoptions.npDA;    % total numbr of grid points
    dx = XmaxDA/npdax; % grid stepsize in x [m]
    dy = YmaxDA/npday; % grid stepsize in y [m]
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

%% Generates or reads lattices with errors
tstart= tic;
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s calcDAdist: starting DA distribution calculation, mode = %s \n', ...
            datetime, mode );
end

if (isempty(fields(ERlat)))
    ERlat = generate_errlatt(RING,ErrorModel,'tunfams',tunfams, ...
    'nseeds',nseeds,'nittune', nittune, 'TolTune', TolTune,...
    'frac', frac, 'useORM0', useORM0f, 'verbose', verboselevel-1);
else
    if (verboselevel>0)
        fprintf('%s calcDAdist: using previously calculated corrected lattices \n', datetime);
    end
    nseeds    = ERlat.inputs.nseeds;
    corrorbf  = ERlat.inputs.corrorbf;
    corrtunf  = ERlat.inputs.corrtunf;
    useORM0f  = ERlat.inputs.useORM0f;
    tunfams   = ERlat.inputs.tunfams;
    nittune   = ERlat.inputs.nittune;
    TolTune   = ERlat.inputs.TolTune; 
    frac      = ERlat.inputs.frac;
end

RINGe        = ERlat.outputs.RINGe;
orb0_stds    = ERlat.outputs.orb0_stds;
orb_stds     = ERlat.outputs.orb_stds;
stab         = ERlat.outputs.stab;
survivalrate = ERlat.outputs.survivalrate;

DAVs= zeros(nang+1,2*(nseeds+1));
DAs = zeros(1,nseeds+1);
DAxdpsp = zeros(npd, nseeds+1);
DAxdpsm = zeros(npd, nseeds+1);
DAydps  = zeros(npd, nseeds+1);
DAVdps  = zeros(3, 2);

%% Checks z0 and nturns for nan
if (not(isempty(ERlat.outputs.rparae{1})))
   rpara = ERlat.outputs.rparae{1};
   etax  = rpara.etax;
   if (isnan(z0))
       if (check_6d(RING))
        z0 = PhysConst.c*(rpara.syncphase-pi)/(2*pi*rpara.revFreq*rpara.harmon); %if z0 not given choose the synchronous phase
       else
        z0=0.0;
       end
       DAoptions.z0=z0;
   end
   %
   % if input number of turns is nan, and lattice is 6d set it to
   %  nsyncT synchrotron period
   if (isnan(nturns))
       if (check_6d(RING))
            nturns = round(nsyncT/rpara.synctune);
       else
            nturns = 1024;
       end
       DAoptions.nturns=nturns;
   end

else
     fprintf('%s Error in calcDAdist: unperturned ring is unstable \n', datetime);
     DAav=nan;
     DAstd=nan;
     DAdist.inputs.RING=RING;
     DAdist.inputs.ErrorModel=ErrorModel;
     DAdist.outputs.DAav=DAav;
     DAdist.outputs.DAstd=DAav;
     DAdist.outputs.DAstd=DAstd;
     DAdist.outputs.orb0_stds=orb0_stds;
     DAdist.outputs.orb_stds=orb_stds;
     DAdist.outputs.RINGe=RINGe;
     DAdist.outputs.rparae=rparae;
     DAdist.outputs.Itunese=Itunesse;
     telapsed=toc(tstart);
     DAdist.outputs.telapsed=telapsed;
     return
end

%% Calculates DAs

if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s Starting DA calculations \n', datetime);
end

for i=1:nseeds+1
 if (verboselevel>0)
    fprintf('%s seed n. %3d \n', datetime, i-1);
 end
 if (stab(i))
  try
   switch mode
       case 'xy'
            if (strcmpi(DAmode,'border')||strcmpi(DAmode,'smart_in'))
                [DAs(i), DAVs(:,2*i-1:2*i)] = ...
                    calcDA_raw(RINGe{i},DAoptions,etax,...
                    rpara.beta0(1),rpara.beta0(2));
            else
                [DAs(i),~] = ...
                    calcDA_raw(RINGe{i},DAoptions,etax,...
                    rpara.beta0(1),rpara.beta0(2));
            end
            if (verboselevel>0)
                fprintf('%s DA = %5.2f mm**2 \n',datetime , DAs(i))
            end

       case 'xydp'
           for j=1:npd
               if (verboselevel>0)
                   fprintf('%s dp = %5.2f %% \n',datetime, dps(j)*100)
               end
               DAoptions.dp=dps(j);
               [~, DAVdps(:,1:2)] = ...
                    calcDA_raw(RINGe{i},DAoptions,etax,...
                    rpara.beta0(1),rpara.beta0(2));
               DAxdpsp(j,i)=DAVdps(1,1);
               DAxdpsm(j,i)=DAVdps(3,1);
               DAydps(j,i)=DAVdps(2,2); 
           end
           DAoptions.dp=dp; 
       otherwise
           fprintf('%s Error in calcDAdist; unknown mode %s \n', datetime, mode);
   end

  catch ME
     fprintf('%s Error in calcDAdist for seed n. %3d \n', datetime, i);
     fprintf('Error message was:%s \n',ME.message);
     DAs(i)=NaN;
  end
 else
     DAs(i)=NaN;
 end
end
DAav     = mean(DAs(2:nseeds+1),'omitnan');
DAstd    = std(DAs(2:nseeds+1),'omitnan');
DAxdppav = mean(DAxdpsp(:,2:nseeds+1),2,'omitnan');
DAxdppst = std(DAxdpsp(:,2:nseeds+1),0,2,'omitnan');
DAxdpmav = mean(DAxdpsm(:,2:nseeds+1),2,'omitnan');
DAxdpmst = std(DAxdpsm(:,2:nseeds+1),0,2,'omitnan');
DAydpav  = mean(DAydps(:,2:nseeds+1),2,'omitnan');
DAydpst  = std(DAydps(:,2:nseeds+1),0,2,'omitnan');

telapsed=toc(tstart);
%% Collects output structure data
DAdist.inputs.RING=RING;
DAdist.inputs.ErrorModel=ErrorModel;
DAdist.inputs.mode=mode;
DAdist.inputs.nseeds=nseeds;
DAdist.inputs.corrorb=corrorbf;
DAdist.inputs.corrtun=corrtunf;
DAdist.inputs.nittune=nittune;
DAdist.inputs.TolTune=TolTune;
DAdist.inputs.tunfams=tunfams;
DAdist.inputs.frac=frac;
DAdist.inputs.useORM0f=useORM0f;

DAdist.outputs.desc=strcat(sprintf('%s',datetime),' : ', desc);
DAdist.outputs.DAoptions=DAoptions;
DAdist.outputs.DAs=DAs;
DAdist.outputs.DAVs=DAVs;
DAdist.outputs.DAav=DAav;
DAdist.outputs.DAstd=DAstd;
DAdist.outputs.orb0_stds=orb0_stds;
DAdist.outputs.orb_stds=orb_stds;
if (fulloutputf)
    DAdist.outputs.RINGe=RINGe;
    DAdist.outputs.rparae=rparae;
    DAdist.outputs.Itunese=Itunese;
    DAdist.outputs.Ftunese=Ftunese;
else
    DAdist.outputs.RINGe={};
    DAdist.outputs.rparae={};
    DAdist.outputs.Itunese={};
    DAdist.outputs.Ftunese={};
end
DAdist.outputs.stab         = stab;
DAdist.outputs.survivalrate = survivalrate;
DAdist.outputs.DAxdpsp  = DAxdpsp;
DAdist.outputs.DAxdppav = DAxdppav;
DAdist.outputs.DAxdppst = DAxdppst;
DAdist.outputs.DAxdpsm  = DAxdpsm;
DAdist.outputs.DAxdpmav = DAxdpmav;
DAdist.outputs.DAxdpmst = DAxdpmst;
DAdist.outputs.DAydps   = DAydps;
DAdist.outputs.DAydpav  = DAydpav;
DAdist.outputs.DAydpst  = DAydpst;

DAdist.outputs.dps  = dps;
DAdist.outputs.telapsed=telapsed;

if(verboselevel>0)
    fprintf('%s DAdist calculation complete \n', datetime);
end

%% Plots DA Distribution and rms orbits
if (plotf)
    if (verboselevel>0)
        fprintf('Plotting DA... \n');
        if (plotorbrmsf)
            plotDAdist(DAdist,'verbose',verboselevel-1,'plotorbrms');
        else
            plotDAdist(DAdist,'verbose',verboselevel-1);
        end
    else
        if (plotorbrmsf)
            plotDAdist(DAdist,'plotorbrms');
        else
            plotDAdist(DAdist);
        end
    end
end
