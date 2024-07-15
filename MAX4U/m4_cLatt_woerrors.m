%% m4_cLatt_woerrors
% This script and the accompanying "m4_cLatt_werrors" run all "cLatt" 
% function options in a series of steps creating/overwriting a structure 
% named "m4T" in the matlab workspace, saving this variable
% in a "m4T.mat" file and saving all matlab output in a logfile
%
% The intended workflow is 
%  1. create a script based on the "m4_cTemplate.m" editing the 
%     statements in the "Lattice Specific Data" section 
%  2. Run the script created in step 1, which calls this script and
%     possibly amd/or also the "m4_cLatt_werrors" script.
%  3. Rename or copy the resulting "m4T" sctructure accordingly.
%
% The "m4_cLatt_werrors" script is a contunuation of this script, 
% containing the more time-consuming Lattice evaluations with errors. 
% Breaking the evaluation into these two parts gives a convenient way 
% of only performing the faster calculations at a first chek of a 
% new lattice.
%% Note 1
% The worskpace variable ACHRO (a cell array with the AT2 lattice 
%         to be evaluated) must be defined before calling this script.
%
%% Note 2 
% The workspace variables below need to be defined prior to 
% running this script, which is typically done by the script created in 
% step 1 above.
%
% lattname                :string containigg the lattice name
% desc                    :descriptive string 
%% Note 3
% The workspace variables below need to be defined prior to running
% this script. If not available certain features of the cLatt function will
% not be available
%
%   cLoptions.ringtune_fams :(1X2) cell array of strings with magnet families
%                           to be uised for tune matching
%   cLoptions.chrom_fams    :(1X2) cell array of strings with magnet families
%                           to be used for chromaticity matching
%   MagnetStrengthLimits : structure with magnet hardware limits
%   ACHROGRD_a1          : cell array with reference AT2 lattice for design
%                          orbit comparison
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
%
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
m4T = cLatt([],'lattname',lattname, 'desc',desc,'cLoptions',cLoptions,...
               'verbose',1);

m4T = cLatt(m4T,'ACHRO',ACHRO,'ACHRO_ref',ACHROGRD_a1,'verbose',1);

m4T = cLatt(m4T,'MagnetStrengthLimits',MagnetStrengthLimits,'verbose',1);

m4T = cLatt(m4T,'basic','verbose',1);

save('m4T', 'm4T');

%
%% Lattice without errors
% Dynamic Aperture
m4T = cLatt(m4T,'DAxy','verbose',2);

save('m4T', 'm4T');


m4T = cLatt(m4T,'DAxydp','verbose',2);

save('m4T', 'm4T');

% Tune Maps
m4T = cLatt(m4T,'TM_xy','TM_gridxy','TM_gridxdp','TM_gridydp','TM_chro','verbose',2);

save('m4T', 'm4T');

m4T = cLatt(m4T,'TM_difxy','TM_difxdp', 'TM_difydp', 'verbose',2);

save('m4T', 'm4T');

% Momentum Aperture
m4T = cLatt(m4T,'LMA','verbose',2);

save('m4T', 'm4T');

%Touschek lifetime
m4T = cLatt(m4T,'TLT','verbose',3);
save('m4T', 'm4T');

