function m4UT=m4Uc_Latt(ACHRO,lattname,desc,cLoptions,ACHRO_ref,MagnetStrengthLimits,varargin)
% Runs all "cLatt" function options in a series of steps 
% creating/overwriting  a structure m4UT, saving this variable
% in a "m4UT.mat" at the end of each step and saving all matlab 
% output in a logfile
%
% The intended workflow is 
%  1. create a script based on the "m4Uc_Template.m" editing the 
%     statements in the "Lattice Specific Data" section 
%
%  2. Run the script created in step 1, which calls this function 
%
% Each time cLatt is called from this function, a progress bar is started 
% with which the calculation may be interrupted for continuation at a 
% later time, adding info to the already existing structure.
% The code below can be used as a template for continuation runs
%
%
%% Inputs
% Mandatory arguments: 
% ACHRO     : AT2 lattice cell array for one achromat. 
% lattname  : string containigg the lattice name
% desc      : descriptive string 
% V0        : total cavity voltage [V], default = 1.8E6; If "auto", then V0
%                                      is calculated to achieve
%                                      a given bukcet height with "VvsBH.m"
% bh        : bucket height, default = "auto". if = "auto", the  
%               bucket height is calculated from V0 using "bheight.m"
% harm      : harmonic number, default = 176)
%
% MagnetStrengthLimits : structure containing hardware limits for magnets 
% ACHRO_ref : reference achromat for calculation of design orbit
%             deviations. if = {} deviations are not calculated
%
% cLoptions : structure containing at least the fields below. Other fields
%             are added at default values. See the full list in the cLatt
%             function.
%   cLoptions.All_famsO : list of all magnet fanilies. If empty, the
%                         function attempst to find all magnet families
%                         in the input lattice                    
%   cLoptions.ringtune_fams :(1X2) cell array of strings with magnet families
%                           to be uised for tune matching
%   cLoptions.chrom_fams    :(1X2) cell array of strings with magnet families
%                           to be used for chromaticity matching
%   cLoptions.eqfam      : (1XN) cell array of strings with 
%                          the names in the MagnetStrengthLmits structure 
%                          that are to be compared for determining 
%                          challenge levels
%   cLoptions.eqsca      :(1XN) array of scaling factors: useful to scale 
%                         the dipole gradients before comparison, 
%                         when the dipole bending angle has been changed 
%                         to compensate for the introduction of 
%                         reverse bends.
% Optional arguments:
% corchro   : if true, chromaticity is corrected
%
% Optional flags
% 'basonly' : if present only lattice structire creation and basic
%             calculations are done


%% History
% PFT 2024/07/06 : first version
% PFT 2024/07/10 : break up of the script into two - without errors and with errors
% PFT 2024/07/13 : turned script into function
% PFT 2024/07/14 : added handling of exitflag from cLatt
% PFT 2024/07/21 : removed calcLMAdist leaving only calcTLTdist (whuc
%                  alrady contains calcLMAdist
%                  incorporated latest change to calcTLT, having fixed
%                  vertical emittance (rather than coupling ratio) 
%                  as a default.
% PFT 2024/07/26 : added calls to generate_errlatt
% PFT 2024/07/27 : updated MAoptions defaults of deltalimit and delta step size
% PFT 2024/07/28 : updated default of DAoptions.DAmode to 'smart_in'                  
% PFT 2024/07/29 : acdded posibility of chromaticity correction
% PFT 2024/07/30 : added input of total Rf voltage
% PFT 2024/07/31 : added basonly flag
%                  added max rms orbt deviation checks of seed survival
%                  improved handling of default values in cLoptions
%                  structure
% PFT 2024/08/07 : changed default value of DAoptions.nturns to nan
% SJ  2024/08/07 : introduced posibility of inputing either voltage or bucket height
%                  to determine RF voltage  and pass this
%                  voltage to RF cavity when generating the RING cell
%                  array
% PFT 2024/08/18 : bug fix , handling of default vaues for DAoptions 

%% Input argument parsing
corchrof         = getoption(varargin,'corchro',false);
V0               = getoption(varargin,'V0',1.8E6);
bh               = getoption(varargin,'bh','auto');
harm             = getoption(varargin,'harm',176);
basonlyf         = any(strcmpi(varargin,'basonly'));

%% General initialisation

if (isfield(cLoptions,'All_fams0'))
    if(isempty(cLoptions.All_famsO))
        Allfams = findFams(ACHRO);
        cLoptions.All_famsO=[Allfams.Dipoles;Allfams.Multipoles];
    end
else
    Allfams = findFams(ACHRO);
    cLoptions.All_famsO=[Allfams.Dipoles;Allfams.Multipoles];
end

%
if (~isfield(cLoptions,'nseeds'))
    cLoptions.nseeds   = 10;     % number of seeds in calculations with errors
end
if (~isfield(cLoptions,'corrorbf')) 
    cLoptions.corrorbf = true;   % if true, corrects the orbit in calculation with errors
end
if (~isfield(cLoptions,'corrtunf')) 
    cLoptions.corrtunf = true;   % if true, corrects the tunes in calculation with errors
end
if (~isfield(cLoptions,'tunfrac')) 
    cLoptions.tunfrac  = 1.0;    % fraction of quad change to be applied aty each step during tune correction
end
if (~isfield(cLoptions,'TolTune')) 
    cLoptions.TolTune  = 1E-3;   % tolerance for tune matching
end
if (~isfield(cLoptions,'nittune')) 
    cLoptions.nittune  = 5;      % max n. of iterations for tune mtaching
end
if (~isfield(cLoptions,'tunfrac')) 
    cLoptions.tunfrac  = 1.0;    % fraction of quad change to be applied aty each step during tune correction
end
if (~isfield(cLoptions,'useORM0f')) 
    cLoptions.useORM0f = true;   % if true, sets the orbit correction to use the orbit respose
%                                  matrix for the unperturbed ring for all iterations
end

if (~isfield(cLoptions,'ErrorModel'))
    cLoptions.ErrorModel = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);
end
%
if (isfield(cLoptions,'DAoptions'))
    if (~isfield(cLoptions.DAoptions,'DAmode'))
        cLoptions.DAoptions.DAmode   = 'smart_in'; % dynamics aperture calculation mode: "border", "grid", "smart_in" or "smart_out"
    end
    if (~isfield(cLoptions.DAoptions,'nturns'))
        cLoptions.DAoptions.nturns   = nan; % number of turns
    end
    if (~isfield(cLoptions.DAoptions,'betax0'))
        cLoptions.DAoptions.betax0   = NaN; % horizontal beta for normalization - if NaN, no normalization is don
    end
    if (~isfield(cLoptions.DAoptions,'betay0'))
        cLoptions.DAoptions.betay0   = NaN; % vertical beta for normalization - if NaN no normalization is done
    end
    if (~isfield(cLoptions.DAoptions,'xmindas'))
        cLoptions.DAoptions.xmindas  = -0.015;% limits of the range in which the DA border is searched in "border" type modes
    end
    if (~isfield(cLoptions.DAoptions,'xmaxdas'))
        cLoptions.DAoptions.xmaxdas  = 0.015;% limits of the range in which the DA border is searched in "border" type modes
    end
    if (~isfield(cLoptions.DAoptions,'ymaxdas'))
        cLoptions.DAoptions.ymaxdas  = 0.007;% limits of the range in which the DA border is searched in "border" type modes
    end
    if (~isfield(cLoptions.DAoptions,'dpmin'))
        cLoptions.DAoptions.dpmin    = -0.04;% minimum dp for xdp/ydp plane caculation
    end
    if (~isfield(cLoptions.DAoptions,'dpmax'))
        cLoptions.DAoptions.dpmax    = 0.04; % maximum dp for xdp/ydp plane caculation
    end
    if (~isfield(cLoptions.DAoptions,'npd'))
        cLoptions.DAoptions.npd      = 11;% number of points along momentum deviation axis
    end
    if (~isfield(cLoptions.DAoptions,'chroms0'))
        cLoptions.DAoptions.chroms0  = [1 1]/20;% Target chromaticity for one superperiod
    end
    if (~isfield(cLoptions.DAoptions,'Tolchrom'))
        cLoptions.DAoptions.TolChrom = [0.0001 0.0001];% Chromaticity tolerances
    end
    if (~isfield(cLoptions.DAoptions,'Nitchro'))
        cLoptions.DAoptions.Nitchro  = 10; % Max n. iterations of chromaticty correction    
    end
    if (~isfield(cLoptions.DAoptions,'dp'))
        cLoptions.DAoptions.dp       = 0.0;% initial dp/p (6d tracking) or fixed dp/p (4d tracking)     
    end
    if (~isfield(cLoptions.DAoptions,'r0'))
        cLoptions.DAoptions.r0      = 0.020; % initial guess for border mode[m]  
    end
    if (~isfield(cLoptions.DAoptions,'nang'))
        cLoptions.DAoptions.nang    = 40;% number of angular steps for border mode      
    end
    if (~isfield(cLoptions.DAoptions,'z0'))
        cLoptions.DAoptions.z0      = nan; % initial longitudinal coordinate (6d tracking). nan uses synchronous phase    
    end
    if (~isfield(cLoptions.DAoptions,'res'))
        cLoptions.DAoptions.res     = 0.0005;  % resolution [m] for border  mode
    end
    if (~isfield(cLoptions.DAoptions,'alpha'))
        cLoptions.DAoptions.alpha   = 1.100;  % da enlargement factor for border search 
    end
    if (~isfield(cLoptions.DAoptions,'XmaxDA'))
        cLoptions.DAoptions.XmaxDA  = 0.015;% Horizontal range is -Xmax to Xmax [m] for "grid" mode
    end
    if (~isfield(cLoptions.DAoptions,'YmaxDA'))
        cLoptions.DAoptions.YmaxDA  = 0.007;% Vertical range is -Xmax to Xmax [m] for "grid" mode
    end
    if (~isfield(cLoptions.DAoptions,'npdax'))
        cLoptions.DAoptions.npdax   = 64; % number of grid points in x direction is 2*npdax+1
    end
    if (~isfield(cLoptions.DAoptions,'npday'))
        cLoptions.DAoptions.npday   = 64; % number of grid points in y direction is  npday+1
    end
else
    cLoptions.DAoptions.DAmode   = 'smart_in';
    cLoptions.DAoptions.nturns   = nan; 
    cLoptions.DAoptions.betax0   = NaN; 
    cLoptions.DAoptions.betay0   = NaN; 
    cLoptions.DAoptions.xmindas  = -0.015;
    cLoptions.DAoptions.xmaxdas  = 0.015;
    cLoptions.DAoptions.ymaxdas  = 0.007;
    cLoptions.DAoptions.dpmin    = -0.04;
    cLoptions.DAoptions.dpmax    = 0.04; 
    cLoptions.DAoptions.npd      = 11;
    cLoptions.DAoptions.chroms0  = [1 1]/20;
    cLoptions.DAoptions.TolChrom = [0.0001 0.0001];
    cLoptions.DAoptions.Nitchro  = 10; 
    cLoptions.DAoptions.dp       = 0.0;
    cLoptions.DAoptions.r0      = 0.020; 
    cLoptions.DAoptions.nang    = 40;
    cLoptions.DAoptions.z0      = nan;
    cLoptions.DAoptions.res     = 0.0005; 
    cLoptions.DAoptions.alpha   = 1.100;  
    cLoptions.DAoptions.XmaxDA  = 0.015;
    cLoptions.DAoptions.YmaxDA  = 0.007;
    cLoptions.DAoptions.npdax   = 64; 
    cLoptions.DAoptions.npday   = 64; 
end
cLoptions.DAoptions.dx = cLoptions.DAoptions.XmaxDA/cLoptions.DAoptions.npdax;  % grid stepsize in x [m]
cLoptions.DAoptions.dy = cLoptions.DAoptions.YmaxDA/cLoptions.DAoptions.npday;  % grid stepsize in y [m]
cLoptions.DAoptions.dxdy = cLoptions.DAoptions.dx*cLoptions.DAoptions.dy; % grid cell area [m**2]
cLoptions.DAoptions.npDA = (2*cLoptions.DAoptions.npdax+1)*(cLoptions.DAoptions.npday+1); %total number of grid points
cLoptions.DAoptions.X0da = zeros(cLoptions.DAoptions.npDA,1);  % horizontal coordinates of grid points [m]
cLoptions.DAoptions.Y0da = zeros(cLoptions.DAoptions.npDA,1);  % vertical coordinates of grid points [m]

k= 1;
for i=0:cLoptions.DAoptions.npday 
    for j=0:2*cLoptions.DAoptions.npdax
        cLoptions.DAoptions.X0da(k) = -cLoptions.DAoptions.XmaxDA+cLoptions.DAoptions.dx*j;
        cLoptions.DAoptions.Y0da(k) =  cLoptions.DAoptions.dy*i;
        k=k+1;
    end
end

%
if (isfield(cLoptions,'TMoptions'))
    if (~isfield(cLoptions.TMoptions,'mode'))
        cLoptions.TMoptions.mode   = 'x';
    end
    if (~isfield(cLoptions.TMoptions,'dp'))
        cLoptions.TMoptions.dp     = 0.0;
    end
    if (~isfield(cLoptions.TMoptions,'npx'))
        cLoptions.TMoptions.npx    = 51;
    end
    if (~isfield(cLoptions.TMoptions,'npy'))
        cLoptions.TMoptions.npy    = 51;
    end
    if (~isfield(cLoptions.TMoptions,'npd'))
        cLoptions.TMoptions.npd    = 51;
    end
    if (~isfield(cLoptions.TMoptions,'xmin'))
        cLoptions.TMoptions.xmin   = -0.010;
    end
    if (~isfield(cLoptions.TMoptions,'xmax'))
        cLoptions.TMoptions.xmax   = +0.010;
    end
    if (~isfield(cLoptions.TMoptions,'ymin'))
        cLoptions.TMoptions.ymin   = 0.0;
    end
    if (~isfield(cLoptions.TMoptions,'ymax'))
        cLoptions.TMoptions.ymax   = 0.007;
    end
    if (~isfield(cLoptions.TMoptions,'dpmin'))
        cLoptions.TMoptions.dpmin  = -0.04;
    end
    if (~isfield(cLoptions.TMoptions,'dpmax'))
        cLoptions.TMoptions.dpmax  = +0.04;
    end
    if (~isfield(cLoptions.TMoptions,'xmin_dm'))
        cLoptions.TMoptions.xmin_dm = -0.012;
    end
    if (~isfield(cLoptions.TMoptions,'xmax_dm'))
        cLoptions.TMoptions.xmax_dm = +0.012;
    end
    if (~isfield(cLoptions.TMoptions,'ymin_dm'))
        cLoptions.TMoptions.ymin_dm = 0.0;
    end
    if (~isfield(cLoptions.TMoptions,'ymax_dm'))
        cLoptions.TMoptions.ymax_dm = 0.007;
    end
    if (~isfield(cLoptions.TMoptions,'dpmin_dm'))
        cLoptions.TMoptions.dpmin_dm = -0.06;
    end
    if (~isfield(cLoptions.TMoptions,'dpmax_dm'))
        cLoptions.TMoptions.dpmax_dm = +0.06;
    end
    if (~isfield(cLoptions.TMoptions,'nturns'))
        cLoptions.TMoptions.nturns   = 1024;
    end
    if (~isfield(cLoptions.TMoptions,'minampx'))
        cLoptions.TMoptions.minampx  = 30E-6;
    end
    if (~isfield(cLoptions.TMoptions,'minampy'))
        cLoptions.TMoptions.minampy  = 30E-6;
    end
    if (~isfield(cLoptions.TMoptions,'method'))
        cLoptions.TMoptions.method   = 4;
    end
    if (~isfield(cLoptions.TMoptions,'smooth'))
        cLoptions.TMoptions.smooth   = false;
    end
else
    cLoptions.TMoptions.mode   = 'x';
    cLoptions.TMoptions.dp     = 0.0;
    cLoptions.TMoptions.npx    = 51;
    cLoptions.TMoptions.npy    = 51;
    cLoptions.TMoptions.npd    = 51;
    cLoptions.TMoptions.xmin   = -0.010;
    cLoptions.TMoptions.xmax   = +0.010;
    cLoptions.TMoptions.ymin   = 0.0;
    cLoptions.TMoptions.ymax   = 0.007;
    cLoptions.TMoptions.dpmin  = -0.04;
    cLoptions.TMoptions.dpmax  = +0.04;
    cLoptions.TMoptions.xmin_dm = -0.012;
    cLoptions.TMoptions.xmax_dm = +0.012;
    cLoptions.TMoptions.ymin_dm = 0.0;
    cLoptions.TMoptions.ymax_dm = 0.007;
    cLoptions.TMoptions.dpmin_dm = -0.06;
    cLoptions.TMoptions.dpmax_dm = +0.06;
    cLoptions.TMoptions.nturns   = 1024;
    cLoptions.TMoptions.minampx  = 30E-6;
    cLoptions.TMoptions.minampy  = 30E-6;
    cLoptions.TMoptions.method   = 4;
    cLoptions.TMoptions.smooth   = false;
end
%
if (isfield(cLoptions,'MAoptions'))
    if (~isfield(cLoptions.MAoptions,'nperiods'))
        cLoptions.MAoptions.nperiods           = 20;
    end
    if (~isfield(cLoptions.MAoptions,'lmafams'))
        cLoptions.MAoptions.lmafams            = 'all';
    end
    if (~isfield(cLoptions.MAoptions,'stepfam'))
        cLoptions.MAoptions.stepfam            = 1;
    end
    if (~isfield(cLoptions.MAoptions,'deltalimit'))
        cLoptions.MAoptions.deltalimit         = 0.3;
    end
    if (~isfield(cLoptions.MAoptions,'initcoord'))
        cLoptions.MAoptions.initcoord          = [0 0 0 0 0 nan]';
    end
    if (~isfield(cLoptions.MAoptions,'delta'))
        cLoptions.MAoptions.delta              = 0.01;
    end
    if (~isfield(cLoptions.MAoptions,'deltastepsize'))
        cLoptions.MAoptions.deltastepsize      = 0.1;
    end
    if (~isfield(cLoptions.MAoptions,'splits'))
        cLoptions.MAoptions.splits             = 10;
    end
    if (~isfield(cLoptions.MAoptions,'split_step_divisor'))
        cLoptions.MAoptions.split_step_divisor = 2;
    end
    if (~isfield(cLoptions.MAoptions,'nturns'))
        cLoptions.MAoptions.nturns             = nan; 
    end
    if (~isfield(cLoptions.MAoptions,'S0min'))
        cLoptions.MAoptions.S0min              = 0.0;
    end
    if (~isfield(cLoptions.MAoptions,'S0max'))
        cLoptions.MAoptions.S0max              = findspos(ACHRO,length(ACHRO)+1);
    end
else
    cLoptions.MAoptions.nperiods           = 20;
    cLoptions.MAoptions.lmafams            = 'all';
    cLoptions.MAoptions.stepfam            = 1;
    cLoptions.MAoptions.deltalimit         = 0.3;
    cLoptions.MAoptions.initcoord          = [1E-5 0 1E-5 0 0 nan]';
    cLoptions.MAoptions.delta              = 0.01;
    cLoptions.MAoptions.deltastepsize      = 0.1;
    cLoptions.MAoptions.splits             = 10;
    cLoptions.MAoptions.split_step_divisor = 2;
    cLoptions.MAoptions.nturns             = nan; 
    cLoptions.MAoptions.S0min              = 0.0;
    cLoptions.MAoptions.S0max              = findspos(ACHRO,length(ACHRO)+1);
end
%
if (isfield(cLoptions,'TLoptions'))
    if(~(isfield(cLoptions.TLoptions,'Ib')))
        cLoptions.TLoptions.Ib                = 0.5/176;
    end
    if(~(isfield(cLoptions.TLoptions,'integrationmethod')))
        cLoptions.TLoptions.integrationmethod = 'integral';
    end
    if(~(isfield(cLoptions.TLoptions,'abstol')))
        cLoptions.TLoptions.abstol            = 1.0e-16;
    end
    if(~(isfield(cLoptions.TLoptions,'reltol')))
        cLoptions.TLoptions.reltol            = 1.0e-16;
    end
    if(~(isfield(cLoptions.TLoptions,'nperiods')))
        cLoptions.TLoptions.nperiods          = 20;
    end
    if(~(isfield(cLoptions.TLoptions,'LMAperiods')))
        cLoptions.TLoptions.LMAperiods        = 1;
    end
    if(~(isfield(cLoptions.TLoptions,'kcoupl')))
        cLoptions.TLoptions.kcoupl            = 'auto';
    end
    if(~(isfield(cLoptions.TLoptions,'emity')))
        cLoptions.TLoptions.emity             = 8.0E-12;
    end
else
    cLoptions.TLoptions.Ib                = 0.5/176;
    cLoptions.TLoptions.integrationmethod = 'integral';
    cLoptions.TLoptions.abstol            = 1.0e-16;
    cLoptions.TLoptions.reltol            = 1.0e-16;
    cLoptions.TLoptions.nperiods          = 20;
    cLoptions.TLoptions.LMAperiods        = 1;
    cLoptions.TLoptions.kcoupl            = 'auto';
    cLoptions.TLoptions.emity             = 8.0E-12;
end

if(isfield(cLoptions,'OCoptions'))
    if(~isfield(cLoptions.OCoptions,'inCOD'))
        cLoptions.OCoptions.inCOD          = [];
    end
    if(~isfield(cLoptions.OCoptions,'neigen'))
        cLoptions.OCoptions.neigen          = [];
    end
    if(~isfield(cLoptions.OCoptions,'cflags'))
        cLoptions.OCoptions.cflags          = [];
    end
    if(~isfield(cLoptions.OCoptions,'scale'))
        cLoptions.OCoptions.scale          = 0.75;
    end
    if(~isfield(cLoptions.OCoptions,'reforbit'))
        cLoptions.OCoptions.reforbit       = [];
    end
    if(~isfield(cLoptions.OCoptions,'steererlimit'))
        cLoptions.OCoptions.steererlimit   = [];
    end
    if(~isfield(cLoptions.OCoptions,'maxrmsx'))
        cLoptions.OCoptions.maxrmsx        = 0.1E-3;
    end
    if(~isfield(cLoptions.OCoptions,'maxrmsy'))
        cLoptions.OCoptions.maxrmsy        = 0.2E-3;
    end
else
    cLoptions.OCoptions.inCOD          = [];
    cLoptions.OCoptions.neigen         = [];
    cLoptions.OCoptions.cflags         = [];
    cLoptions.OCoptions.scale          = 0.75;
    cLoptions.OCoptions.reforbit       = [];
    cLoptions.OCoptions.steererlimit   = [];
    cLoptions.OCoptions.maxrmsx        = 0.1E-3;
    cLoptions.OCoptions.maxrmsy        = 0.2E-3;
end

if(isfield(cLoptions,'GOoptions'))
    if (~isfield(cLoptions.GOoptions,'GOmode'))
        cLoptions.GOoptions.GOmode = 1; 
    end
    if (~isfield(cLoptions.GOoptions,'chamberHAperture'))
        cLoptions.GOoptions.chamberHAperture   = 11.0E-3;
    end
    if (~isfield(cLoptions.GOoptions,'chamberTomagnetGap'))
        cLoptions.GOoptions.chamberTomagnetGap =  0.5E-3; 
    end
    if (~isfield(cLoptions.GOoptions,'chamberThickness'))
        cLoptions.GOoptions.chamberThickness   =  1.0E-3;
    end
    if (~isfield(cLoptions.GOoptions,'chamberShift'))
        cLoptions.GOoptions.chamberShift       =  4.0E-3;
    end
else
    cLoptions.GOoptions.GOmode             = 1; 
    cLoptions.GOoptions.chamberHAperture   = 11.0E-3;
    cLoptions.GOoptions.chamberTomagnetGap =  0.5E-3;
    cLoptions.GOoptions.chamberThickness   =  1.0E-3;
    cLoptions.GOoptions.chamberShift       =  4.0E-3;
end

%% Creates Lattice Structure and runs basic calculations
[m4UT, exitflag] = cLatt([],'lattname',lattname,'desc',desc,'cLoptions',cLoptions,...
                  'ACHRO',ACHRO,'ACHRO_ref',ACHRO_ref,...
                  'MagnetStrengthLimits',MagnetStrengthLimits,...
                  'verbose',1,'corchro', corchrof,'V0',V0,...
                  'bh',bh,'harm',harm,'basic');

%[m4UT, ~] = cLatt(m4UT,'ACHRO',ACHRO,m4U'ACHRO_ref',ACHRO_ref,'verbose',1);

%[m4UT, ~] = cLatt(m4UT,'MagnetStrengthLimits',MagnetStrengthLimits,'verbose',1);

%[m4UT, exitflag] = cLatt(m4UT,'basic','verbose',1);


save('m4UT', 'm4UT');

if (basonlyf)
    return
end

if (strcmpi(exitflag,'cancelled'))
    return
end
%
%% Lattice without errors
% Dynamic Aperture
[m4UT, exitflag] = cLatt(m4UT,'DAxy','verbose',2);

save('m4UT', 'm4UT');

if (strcmpi(exitflag,'cancelled'))
    return
end

[m4UT, exitflag] = cLatt(m4UT,'DAxydp','verbose',2);

save('m4UT', 'm4UT');

if (strcmpi(exitflag,'cancelled'))
    return
end
% Tune Maps
[m4UT, exitflag] = cLatt(m4UT,'TM_xy','TM_gridxy','TM_gridxdp','TM_gridydp','TM_chro','verbose',2);

save('m4UT', 'm4UT');

if (strcmpi(exitflag,'cancelled'))
    return
end

[m4UT, exitflag] = cLatt(m4UT,'TM_difxy','TM_difxdp', 'TM_difydp', 'verbose',2);

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

% Momentum Aperture
[m4UT, exitflag] = cLatt(m4UT,'LMA','verbose',2);

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

%Touschek lifetime
[m4UT, exitflag] = cLatt(m4UT,'TLT','verbose',2);

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

%% Lattice with errors
% Generates lattices with errors and corrects them
[m4UT, exitflag] = cLatt(m4UT,'bascor','verbose',2);

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

% Tune Maps
[m4UT, exitflag] = cLatt(m4UT,'TMdist','verbose',3);

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end
% Dynamic Aperture
[m4UT, exitflag] = cLatt(m4UT,'DAdistxy','verbose',2); 

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

[m4UT,exitflag] = cLatt(m4UT,'DAdistxydp','verbose',2);  

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

% Momentum Aperture
[m4UT,exitflag] = cLatt(m4UT,'LMAdist','verbose',3);  

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

% Touschek lifetime
[m4UT,exitflag] = cLatt(m4UT,'TLTdist','verbose',3);  

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end
