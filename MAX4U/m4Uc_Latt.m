function m4UT=m4Uc_Latt(ACHRO,lattname,desc,cLoptions,ACHRO_ref,MagnetStrengthLimits)
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
%   cLoptions.GOoptions.GOmode = 1; % selects how the chamber and magnet aperture 
%                        geometries are calcualte
%    Mode 1: Vacuum chamber follows the design orbit,
%            Reverse bends are built as symmetric quadrupoles and
%            have their apertures increased by the amount
%            necesary to accomodate the chamber aperture.
%                
%    Mode 2: Vacuum chamber follows the design orbit,
%            Reverse bends are built as asymmetric quadrupoles 
%            so that their aperture does not need to be changed.
%
%    Mode 3: Vacuum chamber is moved by a fixed amount in the
%            horizontal plane and magnet apertures are not changed. 
%            This causes a reduction of physical aperture
%            and makes the design orbit go off-centre in various 
%            magnets, For the reverse bends this is desired as in this
%            case we assume they are symmetric quadrupoles (as in case
%            1). For sextupoles/octupoles this will generate feed-down
%            to be evaluated
%
% Mode 4: Vacuum chamber is moved by a fixed amount in the
%            horizontal plane for blocks U2 to U4 and magnet 
%            apertures are not changed. This causes a reduction of 
%            physical aperture
%            and makes the design orbit go off-centre in various 
%            magnets, For the reverse bends this is desired as in this
%            case we assume they are symmetric quadrupoles (as in case
%            1). For sextupoles/octupoles this will generate feed-down
%            to be evaluated
%
%   cLoptions.GOoptions.chamberHAperture   : horizontal vacuum chamber
%                                  half-aperture, default = 0.011 m
%   cLoptions.GOoptions.chamberTomagnetGap : chamber to magnet gap, 
%                                          default = 0.5E.4 m;
%   cLoptions.GOoptions.chamberThickness   : chamber thickness = 1.0E-3 m;
%   cLoptions.GOoptions.chamberShift       : Radial shift to be applied to  
%                             vacuum chamber in the range covering 
%                             VC3 to VC7,i.e. from U1/BPM-01 to U5/BPM-01,
%                             default = 4.0E-3 m; 

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

%% General initialisation

if(isempty(cLoptions.All_famsO))
    Allfams = findFams(ACHRO);
    cLoptions.All_famsO=[Allfams.Dipoles;Allfams.Multipoles];
end
%
cLoptions.nseeds   = 10;     % number of seeds in calculations with errors
cLoptions.corrorbf = true;   % if true, corrects the orbit in calculation with errors
cLoptions.corrtunf = true;   % if true, corrects the tunes in calculation with errors
cLoptions.tunfrac  = 1.0;    % fraction of quad change to be applied aty each step during tune correction
cLoptions.TolTune  = 1E-3;   % tolerance for tune matching
cLoptions.nittune  = 5;      % max n. of iterations for tune mtaching
cLoptions.tunfrac  = 1.0;    % fraction of quad change to be applied aty each step during tune correction
cLoptions.useORM0f = true;   % if true, sets the orbit correction to use the orbit respose
%                             matrix for the unperturbed ring for all iterations
%                             
cLoptions.DAoptions.DAmode   = 'border'; % dynamics aperture calculation mode: "border", "grid", "smart_in" or "smart_out"
cLoptions.DAoptions.nturns   = 1024; % number of turns
cLoptions.DAoptions.betax0   = NaN; % horizontal beta for normalization - if NaN, no normalization is don
cLoptions.DAoptions.betay0   = NaN; % vertical beta for normalization - if NaN no normalization is done
cLoptions.DAoptions.xmindas  = -0.015;% limits of the range in which the DA border is searched in "border" type modes
cLoptions.DAoptions.xmaxdas  = 0.015;% limits of the range in which the DA border is searched in "border" type modes
cLoptions.DAoptions.ymaxdas  = 0.007;% limits of the range in which the DA border is searched in "border" type modes
cLoptions.DAoptions.dpmin    = -0.04;% minimum dp for xdp/ydp plane caculation
cLoptions.DAoptions.dpmax    = 0.04; % maximum dp for xdp/ydp plane caculation
cLoptions.DAoptions.npd      = 11;% number of points along momentum deviation axis
cLoptions.DAoptions.chroms0  = [1 1]/20;% Target chromaticity for one superperiod
cLoptions.DAoptions.TolChrom = [0.0001 0.0001];% Chromaticity tolerances
cLoptions.DAoptions.Nitchro  = 10; % Max n. iterations of chromaticty correction
cLoptions.DAoptions.dp       = 0.0;% initial dp/p (6d tracking) or fixed dp/p (4d tracking)     
cLoptions.DAoptions.r0      = 0.020; % initial guess for border mode[m]  
cLoptions.DAoptions.nang    = 40;% number of angular steps for border mode      
cLoptions.DAoptions.z0      = nan; % initial longitudinal coordinate (6d tracking). nan uses synchronous phase    
cLoptions.DAoptions.res     = 0.0005;  % resolution [m] for border  mode
cLoptions.DAoptions.alpha   = 1.100;  % da enlargement factor for border search 
cLoptions.DAoptions.XmaxDA  = 0.015;% Horizontal range is -Xmax to Xmax [m] for "grid" mode
cLoptions.DAoptions.YmaxDA  = 0.007;% Vertical range is -Xmax to Xmax [m] for "grid" mode
cLoptions.DAoptions.npdax   = 64; % number of grid points in x direction is 2*npdax+1
cLoptions.DAoptions.npday   = 64; % number of grid points in y direction is  npday+1
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
%
cLoptions.MAoptions.nperiods           = 20;
cLoptions.MAoptions.lmafams            = 'all';
cLoptions.MAoptions.stepfam            = 1;
cLoptions.MAoptions.deltalimit         = 0.3;
cLoptions.MAoptions.initcoord          = [0 0];
cLoptions.MAoptions.delta              = 0.01;
cLoptions.MAoptions.deltastepsize      = 0.001;
cLoptions.MAoptions.splits             = 10;
cLoptions.MAoptions.split_step_divisor = 2;
cLoptions.MAoptions.nturns             = 1024; 
cLoptions.MAoptions.S0min              = 0.0;
cLoptions.MAoptions.S0max              = findspos(ACHRO,length(ACHRO)+1);

%
cLoptions.TLoptions.Ib                = 0.5/176;
cLoptions.TLoptions.integrationmethod = 'integral';
cLoptions.TLoptions.abstol            = 1.0e-16;
cLoptions.TLoptions.reltol            = 1.0e-16;
cLoptions.TLoptions.nperiods          = 20;
cLoptions.TLoptions.LMAperiods        = 1;
cLoptions.TLoptions.kcoupl            = 'auto';
cLoptions.TLoptions.emity             = 8.0E-12;

%% Creates Lattice Structure and runs basic calculations
[m4UT, ~] = cLatt([],'lattname',lattname, 'desc',desc,'cLoptions',cLoptions,...
                  'ACHRO',ACHRO,'ACHRO_ref',ACHRO_ref,...
                  'MagnetStrengthLimits',MagnetStrengthLimits,...
                  'verbose',1);

%[m4UT, ~] = cLatt(m4UT,'ACHRO',ACHRO,m4U'ACHRO_ref',ACHRO_ref,'verbose',1);

%[m4UT, ~] = cLatt(m4UT,'MagnetStrengthLimits',MagnetStrengthLimits,'verbose',1);

[m4UT, exitflag] = cLatt(m4UT,'basic','verbose',1);

save('m4UT', 'm4UT');

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
