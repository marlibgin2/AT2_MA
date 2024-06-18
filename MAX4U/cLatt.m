function LattStruct = cLatt(varargin)
% Generates or updates a structure containing lattice description and 
% lattice performance data. This function can be used as a high-level wrapper
% to the calcDA/calcDAdist, calcLMA/calLMAdist and calcTunemap packages 
% and its output may be used as input to the lattice2LaHF function.
%
%% Inputs
% Optional arguments
% LattSt  : Structure to which the data will be added. If empty, the
%           output structure is created from scratch.
%           Useful to allow adding lattice evaluation results
%           to an already existing structure keeping what has already
%           been calculated, default = []
%
% ACHRO    : AT2 lattice cell array for one achromat. This may be a
%            structure following LatticeOptData.ACHRO (octupoles not included) 
%            or following LatticeOptData.ACHROGRD (octupoles included). If
%            not given, structure must exist in LattSt.
%
% ACHRO_ref: reference lattice cell array for one achromat. This is
%           typically the MAX IV 3 GeV RING lattice. If the input is given
%           or is avaialble in LattSt structure,
%           the deviation of the closed orbit on ACHO wrt to ACHRO_ref is
%           calculated and added to te output structure
%
% RING     : AT2 lattice cell array for the full ring. If not given and
%            LatticeOptData is available the structure is created by the
%            function using the achromat2ring function
%
% lattname : latttice name, default = '';
% desc     : lattice description, default = '';
% split    : factor by which to split elements when calculating the design
%            orbit, default = 1
% V0       : total cavity volatge [V], defaut = 1.8E6
%
% LatticeOptData : structure with optimization data, default = struct;
% LatticeOptData is typically created by m4U.m and contains,
% among others, the following fields: (if not given, defaults are set)
%
%   LatticeOptData.All_famsO  : cell array of string with the names of 
%                               all families including octupoles
%
%   LatticeOptData.eqfam      : cell array of strings with the names of 
%                               magnet families listed in the structure 
%                               MagnetStrengthLimits corresponding to 
%                               the All_famsO for checking 
%                               "challenge levels" with the 
%                               "chalevel.m" function
%
%   LatticeOptData.eqscal      : (1xN) array of scaling factors,useful to 
%                                scale the dipole gradients before
%                                comparison to MagnetStrengthLimits, 
%                                when the dipole bending angle has been 
%                                changed (wrt 3 degrees for unit cell 
%                                dipoles and 1.5 degrees for matching cell
%                                dipoles) to compensate for the 
%                                introduction of reverse bends. The scaling
%                                factor is (dipole bending
%                                angle)/(referencece bend angle - either 3
%                                or 1.5 degrees)
%
%   LatticeOptData.nseeds      : # of seeds for DA calculation with errors.
%   LatticeOptData.chrom_fams  : (1x2) cell array of strings with names
%                                of magnet families to use for 
%                                chromaticity correction
%   LatticeOptData.ringtune_fams: (1x2) cell array of strings with names
%                                of magnet families to use for 
%                                ring tune correction
%
%   LatticeOptData.DAoptions   : structure with DA aperture calculation
%                                options, with fields (amongs others)
%   
%   LatticeOptData.DAoptions.DAmode   :'grid','border' or 'smart' 
%   LatticeOptData.DAoptions.nturns   : number of turns
%   LatticeOptData.DAoptions.xmaxdas  : limits of the range in which the DA
%                                       border is searched [m]
%   LatticeOptData.DAoptions.xmindas  : limits of the range in which the DA
%                                       border is searched [m] 
%   LatticeOptData.DAoptions.ymaxdas  : limits of the range in which the DA
%                                       border is searched [m]
%   LatticeOptData.DAoptions.dpmin    : minimum dp for xdp/ydp plane caculation
%   LatticeOptData.DAoptions.dpmax    : maximum dp for xdp/ydp plane caculation
%   LatticeOptData.DAoptions.npd      : number of points along momentum deviation axis
%   LatticeOptData.DAoptions.TolChrom : (1x2) array of Chromaticity tolerances
%   LatticeOptData.DAoptions.Nitchro  : Max n. iterations of chromaticity correction
%
%   LatticeOptData.TMoptions   : structure with Tune Map calculation
%                                options, with fields (amongs others)
%   LatticeOptData.TMoptions.npx   : number of points along x axis       
%   LatticeOptData.TMoptions.npy   : number of oints along y axis
%   LatticeOptData.TMoptions.npd   : numer of points along momentum
%                                    deviation axis
%   LatticeOptData.TMoptions.xmin  : min horizontal amplitude [m]
%   LatticeOptData.TMoptions.xmax  : max horizontal amplitude [m]
%   LatticeOptData.TMoptions.ymax  : max vertical amplitude [m]
%   LatticeOptData.TMoptions.ymin  : min vertical amplitude [m]
%   LatticeOptData.TMoptions.dpmax : max energy deviation for chromatic tune footprint
%   LatticeOptData.TMoptions.dpmin : min energy deviation for chromatic tune footprint
%   LatticeOptData.TMoptions.nturns: number of turns for tune calculation
%   LatticeOptData.TMoptions.minampx : minimum absolute value of amplitude in horizontal direction,
%   LatticeOptData.TMoptions.minampy : minimum absolute value of amplitude in horizontal direction,
%   LatticeOptData.TMoptions.method  :1: Highest peak in fft
%                                     2: Interpolation on fft results
%                                     3: Windowing + interpolation (default)
%                                     4: NAFF (this is always the method in case the mode is "diff" 
%   LatticeOptData.TMoptions.smooth : if true, uses nearby calcualted tunes
%                                     to avoid unphysical jumps
%
%   LatticeOptData.MAoptions.lmafams: cell array of strings with names of magnet families at which LMA
%                     is to be calculated. If = 'all' then all non-zero length elements are included           = 'all';
%   LatticeOptData.MAoptions.stepfam: specifies only one every stepfam elements are included
%   LatticeOptData.MAoptions.deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
%   LatticeOptData.MAoptions.initcoord: initial coordinates [x0 y0]
%   LatticeOptData.MAoptions.delta: initial guess for momentum aperture 
%   LatticeOptData.MAoptions.deltastepsize: step size for LMA search
%   LatticeOptData.MAoptions.splits: number of iterations of step division
%   LatticeOptData.MAoptions.split_step_divisor: factor to reduce step size at each iteration
%   LatticeOptData.MAoptions.nturns: number of turns. If nan then number of turns is chosen as 1.2/Qs                       = 500; % nan mean use 1.2/Qs
%   LatticeOptData.MAoptions.S0max: maximum longitudinal position at which to calculate LMA             = findspos(LatticeOptData.RING,length(LatticeOptData.RING)+1)/20;
%   LatticeOptData.MAoptions.S0min: minimum longitudinal position at which to calculate LMA             = 0.0;
%
% MagnetStrengthLimits :structure with magnet strength Limits, defaut =
%                                                                struct
%
% verbose              :defines level of verbose output, default=0, i.e. no output 
%
% Optional flags
%
% 'basic'     : calculates basic lattice performance data - atsummary
% 'all'       : performs all calculations
% 'DAxy'      : calculates dynamic aperture on the (x,y) plane on-energy 
%               and off-energy without errors.
% 'DAxydp'    : calculates dynamic aperture on the (x,dp) and (y,dp) planes  
%               without errors.
% 'DAdistxy'  : calculates DA distribution with errors on the (x,y) plane. 
% 'DAdistxydp': calculates DA distribution with errors on the (x,dp) and (y,dp) planes.
% 'TM_xy'     : calculates tune map along x and y axes
% 'TM_gridxy' : calculates tunes on a grid of points on on (x,y) plane.
% 'TM_gridxdp': tunes on a grid of points on the (x,dp) plane
% 'TM_gridydp': calculates tunes on a grid of points on the (x,dp) plane
% 'TM_difxy'  : calculates tune diffusion map on the (x,y) plane 
% 'TM_difxdp' : calculates tune diffusion map on the (x,dp) plane
% 'TM_difydp' : calculates tune diffusion map on the (y,dp) plane 
% 'TM_chro'   : calculates chromatic tune map
% 'LMA'       : calculates Local Momentum Aperture for one achromat wihtout
%               errors
% 'LMAdist'   : calculates Local Momentum Aperture for whole ring with errors 


%% Outputs
% LattStruct is a structure with fields
% LattStruct.Lattice_Name: string with lattice name following naming convention
% LattStruct.Description : string with verbose description of how latice was developed
%
% LattStruct.ACHROMAT    : AT2 lattice cell array for one achromat (an echo of the input)
% LattStruct.ACHROMAT_ref: AT2 lattice cell array for one reference achromat (an echo of the input)
%
% LattStruct.V0          : total Rf voltage (echo of input)
% LattStruc.LattData.LatticeOptData : echo of input.
%
% LattStruc.LattData.XAllO : strengths of all families (includes octupoles)
%
% LattStruc.LattData.RINGGRD  : Full ring latice 6d enabled with same magnet strengths as
%            ACHROMAT. Octupoles are the same as in ACHRO (if they exist
%            there). Otehrwise they are zero.
%
% LattStruct.LattData.ACHROMAT_ZC : achromat for whih the chromatoicity is
%                                   set to zero
% LattStruct.LattData.DesignOrbit: structure with the following fields 
%                          describing the design orbit for one achromat.
% LattStruct.LattData.DesignOrbit.s2d: orbit length
% LattStruct.LattData.DesignOrbit.x2d: horizontal coordinate [m]
% LattStruct.LattData.DesignOrbit.y2d: vertical coordinate [m]
% LattStruct.LattData.DesignOrbit.a2d: angle [rad]
% LattStruct.LattData.DesignOrbit.baa: change in angle [rad]
% LattStruct.LattData.DesignOrbit.ban: Element family if bending angle is not zero
% LattStruct.LattData.DesignOrbit.s2d_ref: orbit length (reference lattice)
% LattStruct.LattData.DesignOrbit.x2d_ref: horizontal coordinate [m]
%                                           (reference lattice)
% LattStruct.LattData.DesignOrbit.y2d_ref: vertical coordinate [m]
%                                           (reference lattice)
% LattStruct.LattData.DesignOrbit.a2d_ref: angle [rad] reference lattice
%
% LattStruct.LattData.DesignOrbit.Deviation: Distance between central
%                                             orbit of the inout achromat
%                                             and the refeence achromat
%
% LattStruct.LattData.MagnetStrengthLimits : Magnet Strength Limits
%                                            structure
% LattStruct.LattData.CLv : Challenge Levels structure, output from
%                           chalevel.
%
% LattStruct.LattPerf.atsummary : output of atsummary run on RINGGRD
%
% Note: The DAs calculated below are for the full ring withe the
% chromaticity as given in the inut lattice
%
% LattStruct.LattPerf.DA.xy_0   : Dynamic aperture on (x,y) plane for
%                                 on-momentum particles without errors
% LattStruct.LattPerf.DA.xy_p3  : Dynamic aperture on (x,y) plane for +3%
%                                 momentum deviation without errors
% LattStruct.LattPerf.DA.xy_m3  : Dynamic aperture on (x,y) plane for -3%
%                                 momentum deviation without errors
% LattStruct.LattPerf.DA.xydp   : Dynamic aperture on (x,dp) and (y,dp) 
%                                 planes without errors
% LattStruct.LattPerf.DAdist.xy : Dynamic aperture on the (x,y) plane fpor
%                                 on-momentum particles with errors
% LattStruct.LattPerf.DAdist.xydp: Dynamic aperture on the (x,dp) and y,dp)
%                                  planes with errors
% Note: The tune maps below are calculated for the zero chromaticity
%       achromat !
% LattStruct.LattPerf.TM.xy: tune map along x and y axis (ADTS)
% LattStruct.LattPerf.TM.gridxy: tune map on a grid of points in (x,y) plane.
% LattStruct.LattPerf.TM.gridxdp: tune map a grid of points in (x,dp) plane
% LattStruct.LattPerf.TM.gridydp: tune map a grid of points in (y,dp) plane
% LattStruct.LattPerf.TM.difxy  : tune diffusion map in the xy plane 
% LattStruct.LattPerf.TM.difxdp : tune diffusion map in the xdp plane
% LattStruct.LattPerf.TM.difydp : tune diffusion map in the ydp plane 
% LattStruct.LattPerf.TM.chro   : tunes vs momentum deviation
%
% LattStructe.Lattperf.LMA      : local momentum aperture without errors
% LattStructe.Lattperf.LMAdist  : local momentum aperture without errors
%% Usage examples
% First creates the structure, allocating all fields ad sets name and
% description
%
% m4U_240316_b03_01_03_01=cLatt('lattname','m4U_230628_b03_01_01_01','desc','blah blah blah');
%
% Next, adds the AT2 lattice cell array to the already existing Lattice
% structure
% m4U_240331_b03_01_04_01=cLatt('LattSt',m4U_240316_b03_01_03_01,'ACHRO',ACHRO_b3);
%
% Now adds Structure with data for DA calculations and Magnet Strength
% Limit Checks. 
%
% m4U_240331_b03_01_04_01=cLatt('LattSt',m4U_240316_b03_01_03_01,'LatticeOptData', LattiuceOptData, 'verbose', 2);
% m4U_240331_b03_01_04_01=cLatt('LattSt',m4U_240316_b03_01_03_01,'MagnetStrengthLimits', MagnetStrengthLimits);
%
% And now adds results of DA calculation with errors (xy plane)
% m4U_240331_b03_01_04_01=cLatt('LattSt',m4U_240316_b03_01_03_01,'DAdistxy','verbose', 2);
%
% And now adds a reference lattice to allow calcuyaltion of the deviation of
% the design orbit
%
% m4U_240331_b03_01_04_01=cLatt('LattSt',m4U_240316_b03_01_03_01,'ACHRO_ref',ACHRO_a1);
%
%% History
% PFT 2024/05 first version
% PFT 2024/06/03 : added info for magnet challenge level calculation
% PFT 2024/06/04 : updated output structure with new sub-fields
% PFT 2024/06/08 : updated output structre with empty fields for items to
%                  be calculated later
% PFT 2024/06/11 : added calculation of central orbit
% PFT 2024/06/16 : added calculatin of LMA and LMA distribution with errors
% PFT 2024/06/18 : upadated extraction of magnet strengths to a more 
%                  general method built into the chalevel function
%
%% Input argument parsing

lattname             = getoption(varargin,'lattname','');
desc                 = getoption(varargin,'desc','');
ACHRO                = getoption(varargin,'ACHRO',{});
ACHRO_ref            = getoption(varargin,'ACHRO_ref',{});
RING                 = getoption(varargin,'RING',{});
LattSt               = getoption(varargin,'LattSt',struct);
LatticeOptData       = getoption(varargin,'LatticeOptData',struct);
MagnetStrengthLimits = getoption(varargin,'MagnetStrengthLimits',struct);
split                = getoption(varargin,'split',1);
V0                   = getoption(varargin,'V0',1.8E6);

verboselevel         = getoption(varargin,'verbose',0);

basicf      = any(strcmpi(varargin,'basic'));
allf        = any(strcmpi(varargin,'all'));
DAxyf       = any(strcmpi(varargin,'DAxy'));
DAxydpf     = any(strcmpi(varargin,'DAxydp'));
DAdistxyf   = any(strcmpi(varargin,'DAdistxy'));
DAdistxydpf = any(strcmpi(varargin,'DAdistxydp'));
TM_xyf      = any(strcmpi(varargin,'TM_xy'));
TM_gridxyf  = any(strcmpi(varargin,'TM_gridxy'));
TM_gridxdpf = any(strcmpi(varargin,'TM_gridxdp'));
TM_gridydpf = any(strcmpi(varargin,'TM_gridydp'));
TM_difxyf   = any(strcmpi(varargin,'TM_difxy'));
TM_difxdpf  = any(strcmpi(varargin,'TM_difxdp'));
TM_difydpf  = any(strcmpi(varargin,'TM_difydp'));
TM_chrof    = any(strcmpi(varargin,'TM_chro'));
LMAf        = any(strcmpi(varargin,'LMA'));
LMAdistf    = any(strcmpi(varargin,'LMAdist'));


%% Constructs output structure template
fprintf(' ************* \n');
fprintf('%s Starting lattice structure creation/update/evaluation. \n', datetime);

if (isempty(fieldnames(LattSt)))
    answer = questdlg('Output Structure will be overwritten !', ...
	        'Warning !', ...
	        'OK','Cancel','Cancel');

    switch answer
        case 'OK'
            if (verboselevel>0)
                fprintf('%s Creating lattice structure from scratch. \n', datetime);
            end

        case 'Cancel'    
            fprintf('%s Aborting... \n', datetime);
            return
    end
    
    LattStruct.Lattice_Name = '';
    LattStruct.Description  = '';
    LattStruct.ACHROMAT = {};
    %
    LattStruct.LattData.ACHROMAT_ref = {};
    LattStruct.LattData.V0 = V0;
    LattStruct.LattData.LatticeOptData=struct;
    LattStruct.LattData.XAllO=[];
    %
    LattStruct.LattData.MagnetStrengthLimits=struct;
    LattStruct.LattData.CLv=struct;
    LattStruct.LattData.ACHROMAT_ZC={};
    LattStruct.LattData.RINGGRD={};
    LattStruct.LattData.DesignOrbit.s2d=[];
    LattStruct.LattData.DesignOrbit.x2d=[];
    LattStruct.LattData.DesignOrbit.y2d=[];
    LattStruct.LattData.DesignOrbit.a2d=[];
    LattStruct.LattData.DesignOrbit.baa=[];
    LattStruct.LattData.DesignOrbit.ban={};
    LattStruct.LattData.DesignOrbit.s2d_ref=[];
    LattStruct.LattData.DesignOrbit.x2d_ref=[];
    LattStruct.LattData.DesignOrbit.y2d_ref=[];
    LattStruct.LattData.DesignOrbit.a2d_ref=[];
    LattStruct.LattData.DesignOrbit.Deviation=[];
    %
    LattStruct.LattPerf.atsummary = struct;
    LattStruct.LattPerf.DA.xy_0   = struct;
    LattStruct.LattPerf.DA.xy_p3  = struct;
    LattStruct.LattPerf.DA.xy_m3  = struct;
    LattStruct.LattPerf.DA.xydp   = struct;
    LattStruct.LattPerf.DAdist.xy   = struct;
    LattStruct.LattPerf.DAdist.xydp = struct;

    LattStruct.LattPerf.TM.xy      = struct;
    LattStruct.LattPerf.TM.xy      = struct;
    LattStruct.LattPerf.TM.gridxy  = struct;
    LattStruct.LattPerf.TM.gridxdp = struct;
    LattStruct.LattPerf.TM.gridydp = struct;
    LattStruct.LattPerf.TM.difxy   = struct;
    LattStruct.LattPerf.TM.difxdp  = struct;
    LattStruct.LattPerf.TM.difydp  = struct;
    LattStruct.LattPerf.TM.chro    = struct;

    LattStruct.LattPerf.LMA        = struct;
    LattStruct.LattPerf.LMAdist    = struct;

else
    LattStruct=LattSt;
    if (verboselevel>0)
        fprintf('%s Lattice data will be added to an existing structure \n', datetime);
    end
end

%
%% Collects output structure lattice data
%
if (not(isempty(lattname)))
    LattStruct.Lattice_Name = lattname;
end
if (not(isempty(desc)))
    LattStruct.Description = desc;
end

if (not(isempty(ACHRO)))
    LattStruct.ACHROMAT = ACHRO;
end

ACHRO = LattStruct.ACHROMAT;

if (not(isempty(RING)))
    LattStruct.LattData.RINGGRD = RING;
end

RING=LattStruct.LattData.RINGGRD;

if (isempty(ACHRO))
    fprintf('%s cLatt Warning: no achromat structure available. Interrupting...\n', datetime);
    return
end

if (not(isempty(ACHRO_ref)))
    LattStruct.LattData.ACHROMAT_ref = ACHRO_ref;
end

ACHRO_ref     = LattStruct.LattData.ACHROMAT_ref;

if (isempty(ACHRO_ref))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: no reference achromat structure available. \n', datetime);
    end
end

if (split>1)
    ACHRO_SP=splitlat(ACHRO,split);
    [s2d, x2d, y2d, a2d, baa, ban] = Survey2D(ACHRO_SP,9.0*pi/180);
else
    [s2d, x2d, y2d, a2d, baa, ban] = Survey2D(ACHRO,9.0*pi/180);
end
   
LattStruct.LattData.DesignOrbit.s2d=s2d;
LattStruct.LattData.DesignOrbit.x2d=x2d;
LattStruct.LattData.DesignOrbit.y2d=y2d;
LattStruct.LattData.DesignOrbit.a2d=a2d;
LattStruct.LattData.DesignOrbit.baa=baa;
LattStruct.LattData.DesignOrbit.ban=ban;

if (not(isempty(ACHRO_ref)))
    if (split>1)
        ACHRO_ref_SP=splitlat(ACHRO_ref,split);
        [s2d_ref, x2d_ref, y2d_ref, a2d_ref,~,~] = Survey2D(ACHRO_ref_SP,9.0*pi/180);
    else
        [s2d_ref, x2d_ref, y2d_ref, a2d_ref,~,~] = Survey2D(ACHRO_ref,9.0*pi/180);
    end
    LattStruct.LattData.DesignOrbit.s2d_ref = s2d_ref;
    LattStruct.LattData.DesignOrbit.x2d_ref = x2d_ref;
    LattStruct.LattData.DesignOrbit.y2d_ref = y2d_ref;
    LattStruct.LattData.DesignOrbit.a2d_ref = a2d_ref;
    [~, ia, ~] = unique(x2d_ref);
    x2d_refu = x2d_ref(ia);
    y2d_refu = y2d_ref(ia);
    dist  = zeros(length(x2d),1);
    for i = 1:length(x2d)
        yinter=interp1(x2d_refu,y2d_refu,x2d(i));
        dist(i)    = abs(yinter - y2d(i));
    end
else
    dist=nan(length(x2d));
end

LattStruct.LattData.DesignOrbit.Deviation = dist;

if (not(isempty(fieldnames(LatticeOptData))))
    LattStruct.LattData.LatticeOptData=LatticeOptData;
end
LatticeOptData=LattStruct.LattData.LatticeOptData;

if (isempty(fieldnames(LatticeOptData)))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: LattOptData structure not available. Using defaults... \n', datetime);
    end
    
    LatticeOptData.nittune      = 5; % max n. of iterations for tune mtaching
    LatticeOptData.TolTune      = 1E-3;% tolerance for tune matching
    LatticeOptData.tunfrac      = 1.0; % fraction of quad change to be applied aty each step during tune correction
    LatticeOptData.chrom_fams   = {'S2_a1';'S5_a1'}; % chromaticty correction
    LatticeOptData.ringtune_fams= {'Q1_a1';'Q2_a1'};% ring tunes matching

    LatticeOptData.All_famsO = {'Q1_a1';'Q2_a1';...
                             'D2_a1';'D1_a1';...
                             'Q3_a1';'Q4_a1';'S1_a1';...
                             'S2_a1';'S3_a1';'S4_a1';...
                             'S5_a1';'O1_a1';'O2_a1';...
                             'O3_a1'}; 

    LatticeOptData.eqfam = {'Qfend_Qdend';'Qfend_Qdend';...
                           'dip';      'dipm';...
                           'Qf_Qfm';    'Qf_Qfm';     'Sdend';...
                           'Sfm';       'Sd';         'Sfi_Sfo';...
                           'Sfi_Sfo';   'Oxx_Oxy';    'Oxx_Oxy';...
                           'Oyy'};

    LatticeOptData.eqsca = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];
  
    LatticeOptData.DAoptions.DAmode   = 'border';
    LatticeOptData.DAoptions.nturns   = 1024;
    LatticeOptData.DAoptions.xmindas  = -0.015;
    LatticeOptData.DAoptions.xmaxdas  = 0.015;
    LatticeOptData.DAoptions.ymaxdas  = 0.007;
    LatticeOptData.DAoptions.dpmin    = -0.04;
    LatticeOptData.DAoptions.dpmax    = 0.04;
    LatticeOptData.DAoptions.npd      = 11;
    LatticeOptData.DAoptions.TolChrom = [0.0001 0.0001];
    LatticeOptData.DAoptions.Nitchro  = 10;
    LatticeOptData.ErrorModel = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);
    LatticeOptData.nseeds     = 10;
    LatticeOptData.ACHRO      = {};
    LatticeOptData.ACHROGRD   = {};
    LatticeOptData.RING       = {};
    LatticeOptData.RINGGRD    = {};
    
    LattStruct.LattData.LatticeOptData=LatticeOptData;
end

if (not(isempty(fieldnames(MagnetStrengthLimits))))
    LattStruct.LattData.MagnetStrengthLimits=MagnetStrengthLimits;
end
MagnetStrengthLimits=LattStruct.LattData.MagnetStrengthLimits;

if (isempty(fieldnames(MagnetStrengthLimits)))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: MagnetStrengthLimits structure not available. \n', datetime);
    end
end

if not(isempty(fieldnames(LatticeOptData)))
    eqfam = LatticeOptData.eqfam;
    eqsca = LatticeOptData.eqsca;
    All_famsO = LatticeOptData.All_famsO;
    CLv = chalevel(MagnetStrengthLimits,'mode','LS',...
                              'ACHRO',ACHRO,'eqfam', eqfam,'eqsca',...
                              eqsca,'Fams',All_famsO,'verbose',verboselevel-1)';  
    LattStruct.LattData.CLv=CLv;
    XAllO=CLv.outputs.X0;
    LattStruct.LattData.XAllO=XAllO;

    if (isempty(RING))
        if (verboselevel>0)
            fprintf('%s cLatt Warning: Full ring input not available, creating ring from achromat...\n',datetime);
        end
        RING = atenable_6d(achromat2ring(ACHRO));
        cavpts  = find(atgetcells(RING, 'FamName', 'CAVITY'));
        if (isempty(cavpts))
            cavpts  = find(atgetcells(RING, 'FamName', 'CAV'));
        end
        RING = atSetRingProperties(RING,'rf_frequency','nominal',...
              'cavpts',cavpts,...
              'rf_voltage',V0*6);
    else
        V0=nan;
        LattStruct.LattData.V0=V0;
    end
    RINGGRD=calculateGirderMaps(RING);
    LattStruct.LattData.RINGGRD=RING;
    
    TolChrom     = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
    Nitchro      = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticity correction
    chrom_fams   = LatticeOptData.chrom_fams;
    try 
        [ACHRO_zc, ~, ~]=fitchroit(ACHRO, chrom_fams, [0 0], Nitchro, TolChrom); 
    catch ME
        fprintf('%s cLatt Error in fitting zero chromaticity: ', datetime);
        fprintf('Error message was:%s \n',ME.message);
        ACHRO_zc={};
    end
    LattStruct.LattData.ACHROMAT_ZC = ACHRO_zc;
end

if (isempty(RINGGRD))
    if(not(isempty(RING)))
        RINGGRD=RING;
    end
end

if(isempty(RINGGRD))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: full ring cell array not available\n', datetime);
    end
end

%% Calculates atsummary for full ring (or achrom,at of ring not available
if (basicf||allf)
    if (verboselevel>0)
        fprintf('%s AT summary calculation. \n', datetime);
    end
    if (isempty(RINGGRD))
        LattStruct.LattPerf.atsummary = atsummary(LattStruct.LattData.ACHRO);
    else
        LattStruct.LattPerf.atsummary = atsummary(LattStruct.LattData.RINGGRD);
    end
end

%% Evaluates DAs and for full RING
if (not(isempty(RINGGRD)))
%% Dynamic aperture for full ring without errors on (x,y) plane
  if (DAxyf||allf)
    if (verboselevel>0)
      fprintf('%s DA calculation: on-momentum. \n', datetime);
    end
    DAS_0 = calcDA(RINGGRD,LatticeOptData.DAoptions, 'mode', 'xy', 'dp', 0.00, 'verbose', verboselevel-1);
    if (verboselevel>0)
      fprintf('%s DA calculation: +3 %%. \n', datetime);
    end
    DAS_p3 = calcDA(RINGGRD,LatticeOptData.DAoptions,'mode', 'xy', 'dp',+0.03, 'verbose', verboselevel-1);
    if (verboselevel>0)
      fprintf('%s DA calculation: -3 %%. \n', datetime);
    end
    DAS_m3 = calcDA(RINGGRD,LatticeOptData.DAoptions,'mode', 'xy', 'dp',-0.03, 'verbose', verboselevel-1);
    LattStruct.LattPerf.DA.xy_0=DAS_0;
    LattStruct.LattPerf.DA.xy_p3=DAS_p3;   
    LattStruct.LattPerf.DA.xy_m3=DAS_m3;
  end
%% Dynamic aperture for full ring without errors on (x,dp) and (y,dp) planes
  if (DAxydpf||allf)
    if (verboselevel>0)
      fprintf('%s DA calculation: (x,dp) plane \n', datetime);
    end
    DAS = calcDA(RINGGRD,LatticeOptData.DAoptions, 'mode', 'xydp','verbose', verboselevel-1);
    LattStruct.LattPerf.DA.xydp=DAS;
  end
%% Dynamic aperture for full ring with errors on xy plane
  if (DAdistxyf||allf)
    if (verboselevel>0)
      fprintf('%s DA calculation with errors: on-momentum. \n', datetime);
    end
    DAdist = calcDAdist(RINGGRD,LatticeOptData.ErrorModel,LatticeOptData.DAoptions,...
             'tunfams',LatticeOptData.ringtune_fams,'mode','xy',...
             'corrorb','corrtun','frac',LatticeOptData.tunfrac,...
             'nseeds',LatticeOptData.nseeds,'verbose', verboselevel-1);
    LattStruct.LattPerf.DAdist.xy = DAdist;
  end

%% Dynamic aperture for full ring with errors on xdp an ydp planes
  if (DAdistxydpf)
    if (verboselevel>0)
      fprintf('%s DA calculation with errors: off-momentum. \n', datetime);
    end
    DAdist = calcDAdist(RINGGRD,LatticeOptData.ErrorModel,  LatticeOptData.DAoptions,...
             'tunfams',LatticeOptData.ringtune_fams,'mode','xydp',...
             'corrorb','corrtun','frac',LatticeOptData.tunfrac,...
             'nseeds',LatticeOptData.nseeds,'verbose', verboselevel-1);
    LattStruct.LattPerf.DAdist.xydp = DAdist;
  end
else
    if (verboselevel>0)
        fprintf('%s cLatt Warning : RINGRD structure not available. \n', datetime);
    end
end
%% Evaluates tune maps for an achromat
if (not(isempty(ACHRO_zc)))
%% Tune map along x and y axes
  if(TM_xyf||allf)
      if (verboselevel>0)
        fprintf('%s Tune Map xy (ADTS) \n', datetime);
      end
      tunemap = calcTuneMap(ACHRO_zc,'mode','xy', ...
                'npx',LatticeOptData.TMoptions.npx,...
                'npy', LatticeOptData.TMoptions.npy,...
                'xmin',LatticeOptData.TMoptions.xmin,...
                'xmax',LatticeOptData.TMoptions.xmax,...
                'ymin',LatticeOptData.TMoptions.ymin,...
                'ymax',LatticeOptData.TMoptions.ymax,...
                'nturns', LatticeOptData.TMoptions.nturns,...
                'minampx',LatticeOptData.TMoptions.minampx, ...
                'minampy',LatticeOptData.TMoptions.minampy,...
                'method',LatticeOptData.TMoptions.method,...;
                'smooth',LatticeOptData.TMoptions.smooth,...;
                'plottype','xy');
      LattStruct.LattPerf.TM.xy=tunemap;
  end

%% Tune map on a grid of points in (x,y) plane.
  if (TM_gridxyf||allf)
       if (verboselevel>0)
        fprintf('%s Tune Map grid (x,y) (ADTS). \n', datetime);
       end
       tunemap = calcTuneMap(ACHRO_zc,'mode','gridxy',...
                 'npx',LatticeOptData.TMoptions.npx,...
                 'npy', LatticeOptData.TMoptions.npy,...
                 'xmin',LatticeOptData.TMoptions.xmin,...
                 'xmax',LatticeOptData.TMoptions.xmax,...
                 'ymin',LatticeOptData.TMoptions.ymin,...
                 'ymax',LatticeOptData.TMoptions.ymax,...
                 'nturns', LatticeOptData.TMoptions.nturns,...
                 'minampx',LatticeOptData.TMoptions.minampx, ...
                 'minampy',LatticeOptData.TMoptions.minampy,...
                 'method',LatticeOptData.TMoptions.method,...;
                 'smooth',LatticeOptData.TMoptions.smooth,...;
                 'plottype','gridxy');
       LattStruct.LattPerf.TM.gridxy=tunemap;
  end
%% Tune map on a gri    d of points in (x,dp) plane
  if (TM_gridxdpf||allf)
    if (verboselevel>0)
        fprintf('%s Tune Map grid (x,dp). \n', datetime);
    end
    tunemap = calcTuneMap(ACHRO_zc,'mode','gridxdp',...
               'npx',LatticeOptData.TMoptions.npx,...
               'npd', LatticeOptData.TMoptions.npd,...
               'xmin',LatticeOptData.TMoptions.xmin,...
               'xmax',LatticeOptData.TMoptions.xmax,...
               'dpmin',LatticeOptData.TMoptions.dpmin,...
               'dpmax',LatticeOptData.TMoptions.dpmax,...
               'nturns', LatticeOptData.TMoptions.nturns,...
               'minampx',LatticeOptData.TMoptions.minampx, ...
               'minampy',LatticeOptData.TMoptions.minampy,...
               'method',LatticeOptData.TMoptions.method,...;
               'smooth',LatticeOptData.TMoptions.smooth,...;
               'plottype','gridxdp');
    LattStruct.LattPerf.TM.gridxdp=tunemap;
  end

%% Tune map on a grid of points in (y,dp) plane
  if (TM_gridydpf||allf)
    if (verboselevel>0)
        fprintf('%s Tune Map grid (y,dp). \n', datetime);
    end
    tunemap = calcTuneMap(ACHRO_zc,'mode','gridydp',...
               'npy',LatticeOptData.TMoptions.npy,...
               'npd', LatticeOptData.TMoptions.npd,...
               'ymin',LatticeOptData.TMoptions.ymin,...
               'ymax',LatticeOptData.TMoptions.ymax,...
               'dpmin',LatticeOptData.TMoptions.dpmin,...
               'dpmax',LatticeOptData.TMoptions.dpmax,...
               'nturns', LatticeOptData.TMoptions.nturns,...
               'minampx',LatticeOptData.TMoptions.minampx, ...
               'minampy',LatticeOptData.TMoptions.minampy,...
               'method',LatticeOptData.TMoptions.method,...;
               'smooth',LatticeOptData.TMoptions.smooth,...;
               'plottype','gridydp');
    LattStruct.LattPerf.TM.gridydp=tunemap;
  end

%% Tune diffusion map on a grid of points in (x,y) plane
  if (TM_difxyf||allf)
    if (verboselevel>0)
        fprintf('%s Tune diffusion map (x,y). \n', datetime);
    end
    tunemap = calcTuneMap(ACHRO_zc,'mode','difxy',...
               'npx',LatticeOptData.TMoptions.npx,...
               'npy', LatticeOptData.TMoptions.npy,...
               'xmin',LatticeOptData.TMoptions.xmin,...
               'xmax',LatticeOptData.TMoptions.xmax,...
               'ymin',LatticeOptData.TMoptions.ymin,...
               'ymax',LatticeOptData.TMoptions.ymax,...
               'nturns', LatticeOptData.TMoptions.nturns,...
               'minampx',LatticeOptData.TMoptions.minampx, ...
               'minampy',LatticeOptData.TMoptions.minampy,...
               'method',LatticeOptData.TMoptions.method,...;
               'smooth',LatticeOptData.TMoptions.smooth,...;
               'plottype','difxy');
    LattStruct.LattPerf.TM.difxy=tunemap;
  end

%% Tune diffusion map on a grid of points in (x,dp) plane
  if (TM_difxdpf||allf)
    if (verboselevel>0)
        fprintf('%s Tune diffusion map (x,dp). \n', datetime);
    end
    tunemap = calcTuneMap(ACHRO_zc,'mode','difxdp',...
               'npx',LatticeOptData.TMoptions.npx,...
               'npd', LatticeOptData.TMoptions.npd,...
               'xmin',LatticeOptData.TMoptions.xmin,...
               'xmax',LatticeOptData.TMoptions.xmax,...
               'dpmin',LatticeOptData.TMoptions.dpmin,...
               'dpmax',LatticeOptData.TMoptions.dpmax,...
               'nturns', LatticeOptData.TMoptions.nturns,...
               'minampx',LatticeOptData.TMoptions.minampx, ...
               'minampy',LatticeOptData.TMoptions.minampy,...
               'method',LatticeOptData.TMoptions.method,...;
               'smooth',LatticeOptData.TMoptions.smooth,...;
               'plottype','difxdp');
    LattStruct.LattPerf.TM.difxdp=tunemap;
  end

%% Tune diffusion map on a grid of points in (y,dp) plane
  if (TM_difydpf||allf)
    if (verboselevel>0)
        fprintf('%s Tune diffusion map (y,dp). \n', datetime);
    end
    tunemap = calcTuneMap(ACHRO_zc,'mode','difydp',...
               'npy',LatticeOptData.TMoptions.npy,...
               'npd', LatticeOptData.TMoptions.npd,...
               'ymin',LatticeOptData.TMoptions.ymin,...
               'ymax',LatticeOptData.TMoptions.ymax,...
               'dpmin',LatticeOptData.TMoptions.dpmin,...
               'dpmax',LatticeOptData.TMoptions.dpmax,...
               'nturns', LatticeOptData.TMoptions.nturns,...
               'minampx',LatticeOptData.TMoptions.minampx, ...
               'minampy',LatticeOptData.TMoptions.minampy,...
               'method',LatticeOptData.TMoptions.method,...;
               'smooth',LatticeOptData.TMoptions.smooth,...;
               'plottype','difydp');
    LattStruct.LattPerf.TM.difydp=tunemap;
  end

%% Chromatic tune map
  if (TM_chrof||allf)
    if (verboselevel>0)
        fprintf('%s Chromatic tune map. \n', datetime);
    end
    tunemap = calcTuneMap(ACHRO_zc,'mode','chro',...
               'npd', LatticeOptData.TMoptions.npd,...
               'dpmin',LatticeOptData.TMoptions.dpmin,...
               'dpmax',LatticeOptData.TMoptions.dpmax, ...
               'nturns', LatticeOptData.TMoptions.nturns,...
               'minampx',LatticeOptData.TMoptions.minampx, ...
               'minampy',LatticeOptData.TMoptions.minampy,...
               'method',LatticeOptData.TMoptions.method,...;
               'smooth',LatticeOptData.TMoptions.smooth,...;
               'plottype','chro');
    LattStruct.LattPerf.TM.chro=tunemap;
  end
else
   fprintf('%s cLatt Warning: ACHRO_zc structure not available. \n', datetime);
end
%% Evaluates LMA for a single achromat without errors
if (not(isempty(ACHRO)))
    if (LMAf)
        LMA=calcLMA(ACHRO,LatticeOptData.MAoptions,'verbose',verboselevel-1);
        LattStruct.LattPerf.LMA=LMA;
    end
else
    if (verboselevel>0)
        fprintf('%s cLatt Warning: ACHRO structure not available for LMA calculation. \n', datetime);
    end
end
%% Evaluates LMA for a full ring with errors
if (not(isempty(RINGGRD)))
    if (LMAdistf)
        LMAdist=calcLMAdist(RINGGRD,LatticeOptData.ErrorModel,...
             LatticeOptData.MAoptions,'verbose',verboselevel-1,...
             'corrorb','corrtun','tunfams',LatticeOptData.ringtune_fams);
        LattStruct.LattPerf.LMAdist=LMAdist;
    end
else
    if (verboselevel>0)
        fprintf('%s Error: RING structure not available for LMAdist calcualtion. \n', datetime);
    end
end

fprintf('%s Lattice structure creation/update/evaluation completed. \n', datetime);
fprintf(' ************* \n');

%% Auxiliary functions
function ACHRO_SP=splitlat(ACHRO,split)
    splits=1/split*ones(1,split);
    k=1;
%    ACHRO_SP = cell(numel(ACHRO)*split,1);
    for i=1:numel(ACHRO)
        if (ACHRO{i}.Length>0.0)
            element=atdivelem(ACHRO{i},splits);
            for l=1:numel(element)
                ACHRO_SP{k}=element{l};
                k=k+1;
            end
        else
            ACHRO_SP{k}=ACHRO{i};
            k=k+1;
        end
    end
   