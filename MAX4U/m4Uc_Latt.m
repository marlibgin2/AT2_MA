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
% ACHRO_ref : reference achromat for calcualtion of desing orbit
%             deviations. if = {} deviations are not calcualted
%
% cLoptions : structure containing at least the fields below
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
%                         reverse bends
%
%% History
% PFT 2024/07/06 : first version
% PFT 2024/07/10 : break up of the script into two - without errors and with errors
% PFT 2024/07/13 : turned script into function
% PFT 2024/07/14 : added handling of exitflag from cLatt
%% General initialisation

if(isempty(cLoptions.All_famsO))
    Allfams = findFams(ACHRO);
    cLoptions.All_famsO=[Allfams.Dipoles;Allfams.Multipoles];
end
%
cLoptions.nseeds   = 10;
cLoptions.corrorb  = true;
cLoptions.corrtun  = true;
cLoptions.tunfrac  = 1.0;
cLoptions.TolTune  = 1E-3;% tolerance for tune matching
cLoptions.tunfrac  = 1.0; % fraction of quad change to be applied aty each step during tune correction
%
cLoptions.DAoptions.DAmode   = 'border';
cLoptions.DAoptions.nturns   = 1024;
cLoptions.DAoptions.betax0   = NaN; 
cLoptions.DAoptions.betay0   = NaN;
cLoptions.DAoptions.xmindas  = -0.015;
cLoptions.DAoptions.xmaxdas  = 0.015;
cLoptions.DAoptions.ymaxdas  = 0.007;
cLoptions.DAoptions.XmaxDA   = 0.015;
cLoptions.DAoptions.YmaxDA   = 0.007;
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
cLoptions.DAoptions.YmaxDA  = 0.006;   
cLoptions.DAoptions.npdax   = 64; 
cLoptions.DAoptions.npday   = 64; 
cLoptions.DAoptions.dx = cLoptions.DAoptions.XmaxDA/cLoptions.DAoptions.npdax;
cLoptions.DAoptions.dy = cLoptions.DAoptions.YmaxDA/cLoptions.DAoptions.npday;   % grid stepsize in y [m]
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
cLoptions.TMoptions.nturns = 1024;
cLoptions.TMoptions.minampx= 30E-6;
cLoptions.TMoptions.minampy= 30E-6;
cLoptions.TMoptions.method = 4;
cLoptions.TMoptions.smooth = false;
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
cLoptions.TLoptions.kcoupl            = 0.025;

%% Creates Lattice Structure and runs basic calculations
[m4UT, ~] = cLatt([],'lattname',lattname, 'desc',desc,'cLoptions',cLoptions,...
               'verbose',1);

[m4UT, ~] = cLatt(m4UT,'ACHRO',ACHRO,'ACHRO_ref',ACHRO_ref,'verbose',1);

[m4UT, ~] = cLatt(m4UT,'MagnetStrengthLimits',MagnetStrengthLimits,'verbose',1);

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
[m4UT, exitflag] = cLatt(m4UT,'TLT','verbose',3);

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end

%% Lattice with errors
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
[m4UT,exitflag] = cLatt(m4UT,'TLTdist','verbose',4);  

save('m4UT', 'm4UT');
if (strcmpi(exitflag,'cancelled'))
    return
end
