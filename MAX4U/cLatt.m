function [LattStruct, exitflag] = cLatt(varargin)
% Generates or updates a structure containing lattice description and 
% lattice performance data. This function can be used as a high-level wrapper
% to the calcDA/calcDAdist,calcLMA/calLMAdist, calTuneMap and calcTLT packages 
% and its output may be used as input to the lattice2LaHF function.
%
%% Inputs
% Mandatory argument
% LattSt :  structure to which the data will be added. If empty, the
%           output structure is created from scratch.
%           Useful to allow adding lattice evaluation results
%           to an already existing structure keeping what has already
%           been calculated, default = struct();
%
% Optional arguments
%
% ACHRO    : AT2 lattice cell array for one achromat. If
%            not given, structure must exist in LattSt.
%
% ACHRO_ref: reference lattice cell array for one achromat. This is
%           typically the MAX IV 3 GeV RING lattice. If the input is given
%           or is avaialble in LattSt structure,
%           the deviation of the closed orbit on ACHO wrt to ACHRO_ref is
%           calculated and added to te output structure
%
% RING     : AT2 lattice cell array for the full ring. If not given and
%            cLoptions is available, the structure is created by the
%            function using the achromat2ring function
%
% lattname : latttice name, default = '';
% desc     : lattice description, default = '';
% split    : factor by which to split elements when calculating the design
%            orbit, default = 1
% V0       : total cavity voltage [V], defaut = 1.8E6
%
% cLoptions : structure with various optional settings, default = struct;
%               cLoptions may be a copy of the LatticeOptData structure 
%               created by the "m4U.m" function or it maybe created 
%               separately as done in the "m4Tc_Template" and 
%               "m4_cLatt scripts". It needs to contain (at least) the 
%               following fields: (if not given, defaults are set)
%
% cLoptions.All_famsO  : cell array of string with the names of 
%                               all families including octupoles
%
% cLoptions.eqfam      : (NX1) cell array of strings with the names of 
%                               magnet families listed in the structure 
%                               MagnetStrengthLimits corresponding to 
%                               the All_famsO for checking 
%                               "challenge levels" with the 
%                               "chalevel.m" function
%
% cLoptions.eqscal      : (1XN) array of scaling factors,useful to 
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
% cLoptions.GOoptions : structure with the fields below describing hoinse for the
%                vacuum chamber and magnet aperures in the horizontal
%                plane.
%                Note:at the oment these do not affect the apertures or 
%                fields used for tracking !
%
% cLoptions.GOoptions.GOmode = 1; % selects how the chamber and magnet aperture 
%                        geometries are calculated
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
% cLoptions.GOoptions.chamberHAperture   : horizontal vacuum chamber
%                                  half-aperture, default = 0.011 m
% cLoptions.GOoptions.chamberTomagnetGap : chamber to magnet gap, 
%                                          default = 0.5E.4 m;
% cLoptions.GOoptions.chamberThickness   : chamber thickness = 1.0E-3 m;
% cLoptions.GOoptions.chamberShift       : Radial shift to be applied to  
%                             vacuum chamber in the range covering 
%                             VC3 to VC7,i.e. from U1/BPM-01 to U5/BPM-01,
%                             default = 4.0E-3 m; 
%
% cLoptions.nseeds      : # of seeds for DA,LMA and TLT calculation 
%                         with errors, default = 10
% cLoptions.chrom_fams  : (1x2) cell array of strings with names
%                              of magnet families to use for 
%                              chromaticity correction, default = 
% cLoptions.ringtune_fams: (1x2) cell array of strings with names
%                                of magnet families to use for 
%                                ring tune correction, default =
% cLoptions.corrorb: if true,perforam orbit correction for 
%                    ring with errors, default = true
% cLoptions.corrtun: if true,perforam orbit correction for 
%                    ring with errors, default = true
%
% cLoptions.RBfams: (1XN ) cell array of strings with names of reverse bend
%                          families
%
% ********************************************************
% cLoptions.DAoptions   : structure with DA aperture calculation
%                                options, with fields: 
%
% cLoptions.DAoptions.DAmode : dynamics aperture calculation mode (, "grid", "smart_in" or "mart_out") , default = "border"
% cLoptions.DAoptions.nturns : number of turns, default = 1024;
% cLoptions.DAoptions.betax0 : horizontal beta for normalization - if
%                                   NaN, no normalization is done, default = NaN
% cLoptions.DAoptions.betay0 : vertical beta for normalization - if
%                                   NaN no normalization is done, default = NaN
% cLoptions.DAoptions.chroms0: Target chromaticities for one
%                                    superperiod, default = [1 1]/20;
% cLoptions.DAoptions.TolChrom: Chromaticity tolerances, default = [1E-4 1E-4]%                                
% cLoptions.DAoptions.Nitchro  : Max n. iterations of chromaticty
%                                     correction. default = 10
% cLoptions.DAoptions.xmaxdas  : limits of the range in which the DA
%                                     border is searched, defaule=0.015 m
% cLoptions.DAoptions.xmindas  : limits of the range in which the DA
%                                     border is searched, defaule = -0.015 m
% cLoptions.DAoptions.ymaxdas  : limits of the range in which the DA
%                                     border is searched, defaule =0.007 m
% cLoptions.DAoptions.dp       : initial dp/p (6d tracking) or fixed
%                                     dp/p (4d tracking), defaule = 0.0
% cLoptions.DAoptions.dpmin    : minimum dp for xdp/ydp planbe
%                                     calculation, defaule =-0.04
% cLoptions.DAoptions.dpmax    : maximum dp for xdp/ydp plane
%                                     caculation, defaule = +0.04
% cLoptions.DAoptions.npd      : number of points along momentum deviation axis
%
% cLoptions.DAoptions.r0       : initial guess for border DA search, default = 0.02 m
% cLoptions.DAoptions.nang     : number of angular steps in DA border search, default = 40
% cLoptions.DAoptions.z0       : initial longitudinal coordinate (6d tracking). nan uses synchrnous phase
% cLoptions.DAoptions.res      : resolution for DA border search,
%                                     default = 5E-4 m
% cLoptions.DAoptions.alpha    : enlargement factor for DA border
%                                     search, default = = 1.100 
% cLoptions.DAoptions.XmaxDA   : Horizontal range for grid DA mode is
%                                     -Xmax to Xmax [m], default = 0.015 m,
% cLoptions.DAoptions.YmaxDA   : Vertical range is 0 to Ymax [m],
%                                     default = 0.007 m
% cLoptions.DAoptions.npdax    : number of grid points in x direction
%                                     is 2*npdax+1, default = 64
% cLoptions.DAoptions.npday    : number of grid points in y direction
%                                     is  npday+1, defaulet = 64
% cLoptions.DAoptions.dx = cLoptions.DAoptions.XmaxDA/cLoptions.DAoptions.npdax; % grid stepsize in x [m]
% cLoptions.DAoptions.dy = cLoptions.DAoptions.YmaxDA/cLoptions.DAoptions.npday;   % grid stepsize in y [m]
% cLoptions.DAoptions.dxdy = cLoptions.DAoptions.dx*cLoptions.DAoptions.dy; % grid cell area [m**2]
%
% *********************************************************************
% cLoptions.TMoptions   : structure with Tune Map calculation
%                         options, with fields:
%
% cLoptions.TMoptions.npx   : number of points along x axis       
% cLoptions.TMoptions.npy   : number of oints along y axis
% cLoptions.TMoptions.npd   : numer of points along momentum
%                                    deviation axis
% cLoptions.TMoptions.xmin  : min horizontal amplitude [m]
% cLoptions.TMoptions.xmax  : max horizontal amplitude [m]
% cLoptions.TMoptions.ymax  : max vertical amplitude [m]
% cLoptions.TMoptions.ymin  : min vertical amplitude [m]
% cLoptions.TMoptions.dp    : momentum deviation
% cLoptions.TMoptions.dpmax : max energy deviation for chromatic tune footprint
% cLoptions.TMoptions.dpmin : min energy deviation for chromatic tune footprint
% cLoptions.TMoptions.nturns: number of turns for tune calculation
% cLoptions.TMoptions.minampx : minimum absolute value of amplitude in horizontal direction,
% cLoptions.TMoptions.minampy : minimum absolute value of amplitude in horizontal direction,
% cLoptions.TMoptions.method  :1: Highest peak in fft
%                              2: Interpolation on fft results
%                              3: Windowing + interpolation (default)
%                              4: NAFF (this is always the method in case the mode is "diff" 
% cLoptions.TMoptions.smooth : if true, uses nearby calcualted tunes
%                                   to avoid unphysical jumps
%
% ********************************************************************************
% cloptions.MAoptions: structure with local momentm aperture calculation
%                         options, with fields:
%
% cLoptions.MAoptions.lmafams: cell array of strings with names of magnet families at which LMA
%                     is to be calculated. If = 'all' then all non-zero length elements are included           = 'all';
% cLoptions.MAoptions.stepfam: specifies only one every stepfam elements are included
% cLoptions.MAoptions.deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
% cLoptions.MAoptions.initcoord: initial coordinates [x0 y0]
% cLoptions.MAoptions.delta: initial guess for momentum aperture 
% cLoptions.MAoptions.deltastepsize: step size for LMA search
% cLoptions.MAoptions.splits: number of iterations of step division
% cLoptions.MAoptions.split_step_divisor: factor to reduce step size at each iteration
% cLoptions.MAoptions.nturns: number of turns. If nan then number of turns
%                             is chosen as 1.2/Qs, default = 1024
% cLoptions.MAoptions.S0max: maximum longitudinal position at which to calculate LMA [m], 
%                            default if RING is available = findspos(cLoptions.RING,length(cLoptions.RING)+1)/20;
%                            default if RING is not available = 528/20
% cLoptions.MAoptions.S0min: minimum longitudinal position at which to
%                            calculate LMA [m], default = 0.0
%
% *******************************************************
% cLoptions.TLoptions : structure with Touschek lifetime calculation
%                         options, with fields:
%
% cLoptions.TLoptions.Ib : current per bunch [A],. default = 0.5/176
% cLoptions.TLoptions.integrationmethod : see calcTLT_raw, default =
%                                               'integral'
% cLoptions.TLoptions.abstol: absolute tolerance, default  = 1.0e-16
% cLoptions.TLoptions.reltol: relative tolerance, default  = 1.0e-16
% cLoptions.TLoptions.Nperiods: number of periods, default = 20
% cLoptions.TLoptions.LMAperiods: number of periods for which to calculate the LMA, default = 1;
% cLoptions.TLoptions.kcoupl: coupling ratio, if = 'auto', calculated from the
%                                   vertical emittance, default = 'auto'
% cLoptions.TLoptions.emity : vertical emittance, if kcoupl not 'auto', 
%                             calculated from the coupling ratio, 
%                             default =8.0E-12 [m rad]
%
% *******************************************************
%
% MagnetStrengthLimits :structure with magnet strength Limits, defaut =
%                                                               struct
% verbose              :defines level of verbose output, default=0, i.e. no output 
%
%% Optional flags
%
% 'basic'     : calculates basic lattice performance data - atsummary
% 'all'       : performs all calculations
% 'DAs'       : performs dynamic aperture calculations calculations
% 'TMs'       : performs tunemap calculations
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
% 'TLT'       : calculates Touschek lifetime for one achromat without
%               errors
% 'TLTdist'   : calculates Touschek lifetime for ring with errors
%
% 'TMdist'    : Tune map(diffusion in xy plane) for ring with errors
%% Outputs
% LattStruct is a structure with fields
% LattStruct.Lattice_Name: string with lattice name following naming convention
% LattStruct.Description : string with verbose description of how latice was developed
% LattStruct.ACHROMAT    : AT2 lattice cell array for one achromat (an echo of the input)
% LattStuct.Log          : cell array of strings with log of structure construction 
% 
% LattStruc.cLoptions : echo of input cLoptions structure or set of defaults
%
% *******************************************************************
%
% LattStruct.LattData.V0          : total Rf voltage (echo of input)
% LattStruc.LattData.XAllO   :strengths of all families (includes octupoles)
% LattStruct.LattData.ACHROMAT_ref: AT2 lattice cell array for one reference achromat (an echo of the input)
% LattStruc.LattData.RINGGRD  : Full ring latice 6d enabled with same magnet strengths as
%            ACHROMAT. Octupoles are the same as in ACHRO (if they exist
%            there). Otehrwise they are zero.
%
% LattStruct.LattData.ACHROMAT_ZC : achromat for whih the chromatoicity is
%                                   set to zero
%
% LattStruct.LattData.geometry: structure with the following fields 
%                          describing the design orbit, chamber and magnet
%                          horizontal apertures.
%
% LattStruct.LattData.geometry.DesignOrbit.s2d: orbit length
% LattStruct.LattData.geometry.DesignOrbit.x2d: horizontal coordinate [m]
% LattStruct.LattData.geometry.DesignOrbit.y2d: vertical coordinate [m]
% LattStruct.LattData.geometry.DesignOrbit.a2d: angle [rad]
% LattStruct.LattData.geometry.DesignOrbit.baa: change in angle [rad]
% LattStruct.LattData.geometry.DesignOrbit.ban: Element family if bending angle is not zero
%
% LattStruct.LattData.geometry.ref.DesignOrbit.s2d: orbit length (reference lattice)
% LattStruct.LattData.geometry.ref.DesignOrbit.x2d: horizontal coordinate [m]
%                                           (reference lattice)
% LattStruct.LattData.geometry.ref.DesignOrbit.y2d: vertical coordinate [m]
%                                           (reference lattice)
% LattStruct.LattData.geometry.ref.DesignOrbit.a2d_ref: angle [rad] 
%                                           (reference lattice)
%
% LattStruct.LattData.geometry.DesignOrbit.Deviation: distance between central
%                                             orbit of the input achromat
%                                             and the reference achromat
%
% LattStruct.LattData.geometry.ref.Magnets.HAperture : reference horizontal
%                                                      magnet aperture
% LattStruct.LattData.geometry.ref.Chambers.HApertur : reference horizontal
%                                                      chamber aperture
% LattStruct.LattData.geometry.ref.Magnets.x2d_up                   
% LattStruct.LattData.geometry.ref.Magnets.y2d_up
% LattStruct.LattData.geometry.ref.Magnets.x2d_down
% LattStruct.LattData.geometry.ref.Magnets.y2d_down
% LattStruct.LattData.geometry.ref.Chambers.x2d_up
% LattStruct.LattData.geometry.ref.Chambers.y2d_up
% LattStruct.LattData.geometry.ref.Chambers.x2d_down
% LattStruct.LattData.geometry.ref.Chambers.y2d_down
% LattStruct.LattData.geometry.ref.chThick
%
% LattStruct.LattData.geometry.Magnets.Centre.x2d
% LattStruct.LattData.geometry.Magnets.Centre.y2d
% LattStruct.LattData.geometry.Magnets.Centre.Deviation
% LattStruct.LattData.geometry.Magnets.HAperture
% LattStruct.LattData.geometry.Magnets.Walls.x2d_up
% LattStruct.LattData.geometry.Magnets.Walls.y2d_up
% LattStruct.LattData.geometry.Magnets.Walls.x2d_down
% LattStruct.LattData.geometry.Magnets.Walls.y2d_down
% LattStruct.LattData.geometry.Magnets.Walls.dev_up
% LattStruct.LattData.geometry.Magnets.Walls.dev_down
%
% LattStruct.LattData.geometry.Magnets.magceTobeam
%
% LattStruct.LattData.geometry.Chambers.Centre.x2d
% LattStruct.LattData.geometry.Chambers.Centre.y2d
% LattStruct.LattData.geometry.Chambers.Centre.Deviation
% LattStruct.LattData.geometry.Chambers.HAperture
% LattStruct.LattData.geometry.Chambers.Walls.x2d_up
% LattStruct.LattData.geometry.Chambers.Walls.y2d_up
% LattStruct.LattData.geometry.Chambers.Walls.x2d_down 
% LattStruct.LattData.geometry.Chambers.Walls.y2d_down 
% LattStruct.LattData.geometry.Chambers.Walls.dev_up   
% LattStruct.LattData.geometry.Chambers.Walls.dev_down 
% LattStruct.LattData.geometry.Chambers.chaceTobeam 
% LattStruct.LattData.geometry.Chambers.effectiveAperture
%
% LattStruct.LattData.MagnetStrengthLimits : Magnet Strength Limits
%                                            structure
% LattStruct.LattData.CLv : Challenge Levels structure, output from
%                           chalevel.
% LattStruct.LattData.famLayout : cell array with two coluns. The first 
%              contains strings with names of magnet families in
%              the order they appear in the lattice and the 
%              second contains the corresponding multipolar strengths 
%              of tose elentes - the PolynomB array field of the element.
% ******************************************************************
% LattStruct.LattPerf.atsummary : output of atsummary run on RINGGRD
%
% Note: The DAs calculated below are for the full ring with the
%        chromaticity as given in the input lattice. Each output field
%        is itself a structure with several fields (see calcDA for
%        details).
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
%
% ********************************************************************
% Note: The tune maps below are calculated for the zero chromaticity
%       achromat. Each outputfield is itself a structure with 
%       several fields (see calcTuneMap.m for details)
%
% LattStruct.LattPerf.TM.xy: tune map along x and y axis (ADTS)
% LattStruct.LattPerf.TM.gridxy: tune map on a grid of points in (x,y) plane.
% LattStruct.LattPerf.TM.gridxdp: tune map a grid of points in (x,dp) plane
% LattStruct.LattPerf.TM.gridydp: tune map a grid of points in (y,dp) plane
% LattStruct.LattPerf.TM.difxy  : tune diffusion map in the xy plane 
% LattStruct.LattPerf.TM.difxdp : tune diffusion map in the xdp plane
% LattStruct.LattPerf.TM.difydp : tune diffusion map in the ydp plane 
% LattStruct.LattPerf.TM.chro   : tunes vs momentum deviation
%
% ********************************************************************
% Note: The local momentum apperture and Touschek lifetime output fields
%       below are themselves structures with several fields. For details
%       see calcLMA.m, calcTLT.m, calcLMAdist.m and calcTLTdist.m 
%       for details.
%
% LattStructe.Lattperf.LMA      : local momentum aperture without errors
%
% LattStructe.Lattperf.LMAdist  : local momentum aperture without errors
%
% LattStructe.Lattperf.TL       : Touschek lifetime for achromat without
%                                 errors 
%
% LattStructe.Lattperf.TLdist  : (1x2) array : Average and standard
%                                 deviation of Touscke lifeftime for ring
%                                 with errors [hr]. 
% LattStructe.Lattperf.TLdist  : tune diffusion map for lattice wirh errors
%
%% Usage examples
% First creates the structure, allocating all fields and sets name and
% description
%
% m4UT=cLatt([],'lattname','m4U_230628_b03_01_01_01','desc','blah blah blah');
%
% Next, adds the AT2 lattice cell array to the already existing Lattice
% structure
% m4UT=cLatt(m4UT,'ACHRO',ACHRO_b3);
%
% Now adds Structure with data for DA calculations and Magnet Strength
% Limit Checks. 
%
% m4UT=cLatt(m4UT,'cLoptions', LattiuceOptData, 'verbose', 2);
% m4UT=cLatt(m4UT,'MagnetStrengthLimits', MagnetStrengthLimits);
%
% And now adds results of DA calculation with errors (xy plane)
% m4UT=cLatt(m4UT,'DAdistxy','verbose', 2);
%
% And now adds a reference lattice to allow calcuyaltion of the deviation of
% the design orbit
%
% m4UT=cLatt(m4UT,'ACHRO_ref',ACHRO_a1);
%

%% History
% PFT 2024/05 first version
% PFT 2024/06/03 : added info for magnet challenge level calculation
% PFT 2024/06/04 : updated output structure with new sub-fields
% PFT 2024/06/08 : updated output structre with empty fields for items to
%                  be calculated later
% PFT 2024/06/11 : added calculation of central orbit
% PFT 2024/06/16 : added calculatin of LMA and LMA distribution with errors
% PFT 2024/06/18 : updated extraction of magnet strengths to a more 
%                  general method built into the chalevel function
% PFT 2024/06/22 : renamed LatticeOptData as cLoptions and rechecked that 
%                  all required default values are available.
% PFT 2024/06/26 : further restructured the output structure and
%                  added Touschek lifetime for achromat without errors 
%                  and full ring with errors.
% PFT 2024/06/28 : added calculation of fields and gradients.
% PFT 2024/06/30 : changed handling of 'basic' calculation mode
% PFT 2024/07/06 : added choice of calculation subsets (DAs, TMs)
% PFT 2024/07/07 : added checks for unavailable cLoption fields
%                  changed SetBPMWeights to handle lattice without BPMs
% PFT 2024/07/09 : changed SetBPMWeights to handle lattice with mons or
%                  BPMs and print out a warning if none can be found)
%                : added magnet centre calculation
% PFT 2024/07/12 : separate input argumets for tune maps and diffusion maps 
% PFT 2024/07/13 : added calls to TMdist - tune map with errors
%                  added progress bar
% PFT 2024/07/14 : added exitflag output
% PFT 2024/07/17 : added lattice layout determination
% PFT 2024/07/18 : added calcualtion of chamber and magnet aperture
%                  geometries
% PFT 2024/07/21 : added log of structure creation,
%                  changed LMA calculation without errors from single
%                  achromat 4D to full ring 6D.
% PFT 2024/07/22 : added possibilty of fixing vertical emittance instead
%                  of coupling.
%                  added possibility of using existing LMA structure as
%                  input .
%
%% Input argument parsing
LattSt               = getargs(varargin,[]);

lattname             = getoption(varargin,'lattname','');
desc                 = getoption(varargin,'desc','');
ACHRO                = getoption(varargin,'ACHRO',{});
ACHRO_ref            = getoption(varargin,'ACHRO_ref',{});
RING                 = getoption(varargin,'RING',{});
cLoptions            = getoption(varargin,'cLoptions',struct);
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
DAsf        = any(strcmpi(varargin,'DAs'));
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
TLTf        = any(strcmpi(varargin,'TLT'));
TLTdistf    = any(strcmpi(varargin,'TLTdist'));
TMsf        = any(strcmpi(varargin,'TMs'));
TM_distf    = any(strcmpi(varargin,'TMdist'));

%% Constructs output structure template
fprintf(' ************* \n');
fprintf('%s Starting lattice structure creation/update/evaluation. \n', datetime);

if (isempty(LattSt))
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
            exitflag ='abort';
            return
    end
    
    LattStruct.Lattice_Name = lattname;
    LattStruct.Description  = desc;
    LattStruct.ACHROMAT     = ACHRO;
    if (nargin>1)
        if(ischar(varargin{2}))
            args=varargin{2};
        end
    else
        args='';
    end
    for i=3:nargin
        if (ischar(varargin{i}))
            args=[args,',',varargin{i}];
        end
    end
    LattStruct.Log = {strcat(sprintf('%s', datetime),{': Structure Created, args = '}, {args})};
    LattStruct.cLoptions    = cLoptions;
    %
    LattStruct.LattData.ACHROMAT_ref = ACHRO_ref;
    LattStruct.LattData.V0 = V0;
    LattStruct.LattData.XAllO=[];
    %
    LattStruct.LattData.MagnetStrengthLimits=MagnetStrengthLimits;
    LattStruct.LattData.CLv=struct;
    LattStruct.LattData.ACHROMAT_ZC={};
    LattStruct.LattData.RINGGRD=RING;
    %
    LattStruct.LattData.geometry=struct();

    LattStruct.LattData.FG={};
    LattStruct.LattData.famLayout={};

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

    LattStruct.LattPerf.TL         = struct;
    LattStruct.LattPerf.TLdist     = struct;

    LattStruct.LattPerf.TMdist     = struct;

else
    LattStruct=LattSt;
    if (nargin>1)
        if(ischar(varargin{2}))
            args=varargin{2};
        end
    else
        args='';
    end
    for i=3:nargin
        if (ischar(varargin{i}))
            args=[args,',',varargin{i}];
        end
    end
    if (isfield(LattStruct,'Log'))
        LattStruct.Log=[LattStruct.Log; {strcat(sprintf('%s', datetime),{': args = '}, {args} )}];
    else
        LattStruct.Log={strcat(sprintf('%s', datetime),{': args = '}, {args})};
    end
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
lattname = LattStruct.Lattice_Name;

if (not(isempty(desc)))
    LattStruct.Description = desc;
end
desc = LattStruct.Description;

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
    exitflag='abort';
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

if (not(isempty(fieldnames(MagnetStrengthLimits))))
    LattStruct.LattData.MagnetStrengthLimits=MagnetStrengthLimits;
end
MagnetStrengthLimits=LattStruct.LattData.MagnetStrengthLimits;

if (isempty(fieldnames(MagnetStrengthLimits)))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: MagnetStrengthLimits structure not available. \n', datetime);
    end
end

if (not(isempty(fieldnames(cLoptions))))
    LattStruct.cLoptions=cLoptions;
end

if (isfield(LattStruct,'cLoptions'))
    cLoptions=LattStruct.cLoptions;
end


%% Sets default values for cLoptions structure if not available

if (isempty(fieldnames(cLoptions)))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: cLoptions structure not available. Using defaults... \n', datetime);
    end
    
    cLoptions.nittune       = 5; 
    cLoptions.TolTune       = 1E-3;
    cLoptions.tunfrac       = 1.0; % fraction of quad change to be applied aty each step during tune correction
    cLoptions.nseeds        = 10;
    cLoptions.corrorb       = true;
    cLoptions.corrtun       = true;
    cLoptions.ACHRO         = {};
    cLoptions.ACHROGRD      = {};
    cLoptions.RING          = {};
    cLoptions.RINGGRD       = {};
    cLoptions.chrom_fams    = {'S3_b3';'S5_b3'}; % chromaticty correction
    cLoptions.ringtune_fams = {'Q1_b3';'Q2_b3'};% ring tunes matching
    cLoptions.sext_fams     = {'S1_b3';'S2_b3';'S3_b3';'S4_b3';'S5_b3'};

    cLoptions.All_famsO = {'Q1_b3';'Q2_b3';'R1_b3';...
                             'D2_b3';'D3_b3';'D1_b3';...
                             'Q3_b3';'Q4_b3';'S1_b3';...
                             'S2_b3';'S3_b3';'S4_b3';...
                             'S5_b3';'O1_b3';'O2_b3';...
                             'O3_b3'}; 

    cLoptions.eqfam = {'Qfend_Qdend';'Qfend_Qdend';'Qf_Qfm';...
                           'dip';      'dip';        'dipm';...
                           'Qf_Qfm';    'Qf_Qfm';     'Sdend';...
                           'Sfm';       'Sd';         'Sfi_Sfo';...
                           'Sfi_Sfo';   'Oxx_Oxy';    'Oxx_Oxy';...
                           'Oyy'};

    cLoptions.eqsca  = [1 1 1 1 (1.5+2*3.49E-3*180/pi)/1.5 1 1 1 1 1 1 1 1 1 1 1];
    cLoptions.RBfams = {'R1_b3'};

    cLoptions.GOoptions.GOmode   = 1; 
    cLoptions.GOoptions.chamberHAperture   = 11.0E-3;
    cLoptions.GOoptions.chamberTomagnetGap =  0.5E-3;
    cLoptions.GOoptions.chamberThickness   =  1.0E-3;
    cLoptions.GOoptions.chamberShift       =  0.5E-3;

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
    cLoptions.DAoptions.dx      = cLoptions.DAoptions.XmaxDA/cLoptions.DAoptions.npdax;
    cLoptions.DAoptions.dy      = cLoptions.DAoptions.YmaxDA/cLoptions.DAoptions.npday;   % grid stepsize in y [m]
    cLoptions.DAoptions.dxdy    = cLoptions.DAoptions.dx*cLoptions.DAoptions.dy; % grid cell area [m**2]
    
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
    cLoptions.TMoptions.npx = 11;
    cLoptions.TMoptions.npy = 11;
    cLoptions.TMoptions.npd = 11;
    cLoptions.TMoptions.xmin = -0.006;
    cLoptions.TMoptions.xmax = +0.006;
    cLoptions.TMoptions.ymin = 0.0;
    cLoptions.TMoptions.ymax = 0.004;
    cLoptions.TMoptions.dpmin = -0.04;
    cLoptions.TMoptions.dpmax = +0.04;
    cLoptions.TMoptions.xmin_dm = -0.006;
    cLoptions.TMoptions.xmax_dm = +0.006;
    cLoptions.TMoptions.ymin_dm = 0.0;
    cLoptions.TMoptions.ymax_dm = 0.004;
    cLoptions.TMoptions.dpmin_dm = -0.04;
    cLoptions.TMoptions.dpmax_dm = +0.04;
    cLoptions.TMoptions.dp    = 0.0;
    cLoptions.TMoptions.nturns = 1024;
    cLoptions.TMoptions.minampx = 30E-6;
    cLoptions.TMoptions.minampy = 30E-6;
    cLoptions.TMoptions.method = 4;
    cLoptions.TMoptions.smooth = false;
%
    cLoptions.MAoptions.lmafams            = 'all';
    cLoptions.MAoptions.stepfam            = 1;
    cLoptions.MAoptions.deltalimit         = 0.3;
    cLoptions.MAoptions.initcoord          = [0 0];
    cLoptions.MAoptions.delta              = 0.01;
    cLoptions.MAoptions.deltastepsize      = 0.001;
    cLoptions.MAoptions.splits             = 10;
    cLoptions.MAoptions.split_step_divisor = 2;
    cLoptions.MAoptions.nturns             = 1024; 
%
    cLoptions.TLoptions.Ib                = 0.5/176;
    cLoptions.TLoptions.integrationmethod = 'integral';
    cLoptions.TLoptions.abstol            = 1.0e-16;
    cLoptions.TLoptions.reltol            = 1.0e-16;
    cLoptions.TLoptions.Nperiods          = 20;
    cLoptions.TLoptions.LMAperiods        = 1;
    cLoptions.TLoptions.kcoupl            = 'auto';
    cLoptions.TLoptions.emity             = 8.0E-12;

    if (not(isempty(cLoptions.RING)))
        cLoptions.MAoptions.S0max         = findspos(cLoptions.RING,length(cLoptions.RING)+1)/20;
    else
        cLoptions.MAoptions.S0max         = 528/20;
    end

    cLoptions.MAoptions.S0min             = 0.0;

    cLoptions.ErrorModel = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);
    
    LattStruct.cLoptions=cLoptions;
end

%% Wait bar to allow for interruptions without loss of data
fb=waitbar(0,'Lattice Structure Creation/Update', 'Name','Progress', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(fb,'canceling',0);
frac=0.0;
dfrac=100/11;
%% Lattice family layout
if (basicf||allf)
    if (verboselevel>0)
        fprintf('%s cLatt: Determining lattice family layout \n', datetime);
    end
    fLO = famLayout(ACHRO);
    LattStruct.LattData.famLayout=fLO;
end
%% Calculates design orbit
if (basicf||allf)
    if (verboselevel>0)
        fprintf('%s cLatt: calculating design orbit \n', datetime);
    end
    if (split>1)
        ACHRO_SP     = splitlat(ACHRO,split);
    else
        ACHRO_SP     = ACHRO;
    end

    [s2d, x2d, y2d, a2d, baa, ban] = Survey2D(ACHRO_SP,9.0*pi/180);
    
    npo=numel(x2d);

    LattStruct.LattData.geometry.DesignOrbit.s2d=s2d;
    LattStruct.LattData.geometry.DesignOrbit.x2d=x2d;
    LattStruct.LattData.geometry.DesignOrbit.y2d=y2d;
    LattStruct.LattData.geometry.DesignOrbit.a2d=a2d;
    LattStruct.LattData.geometry.DesignOrbit.baa=baa;
    LattStruct.LattData.geometry.DesignOrbit.ban=ban;

    if (not(isempty(ACHRO_ref)))
        if (split>1)
            ACHRO_ref_SP=splitlat(ACHRO_ref,split);
        else
            ACHRO_ref_SP=ACHRO_ref;
        end
        
        [s2d_ref, x2d_ref, y2d_ref, a2d_ref,~,~] = Survey2D(ACHRO_ref_SP,9.0*pi/180);
        npo_ref = numel(x2d_ref);

        LattStruct.LattData.geometry.ref.DesignOrbit.s2d = s2d_ref;
        LattStruct.LattData.geometry.ref.DesignOrbit.x2d = x2d_ref;
        LattStruct.LattData.geometry.ref.DesignOrbit.y2d = y2d_ref;
        LattStruct.LattData.geometry.ref.DesignOrbit.a2d = a2d_ref;
        [~, ia, ~] = unique(x2d_ref);
        x2d_refu = x2d_ref(ia);
        y2d_refu = y2d_ref(ia);
        dist  = zeros(npo,1);
        for i = 1:npo
            yinter=interp1(x2d_refu,y2d_refu,x2d(i));
            dist(i)    = abs(yinter - y2d(i));
        end
    else
        dist=nan(length(x2d),1);
    end

    LattStruct.LattData.geometry.DesignOrbit.Deviation = dist;
end
%% Calculates magnet and chamber geometry
if (basicf||allf)
    
    chamberHAperture_ref   = cLoptions.GOoptions.chamberHAperture;
    chamberTomagnetGap     = cLoptions.GOoptions.chamberTomagnetGap; 
    chamberThickness       = cLoptions.GOoptions.chamberThickness;
    magHAperture_ref       = chamberHAperture_ref+chamberTomagnetGap +...
                             chamberThickness;
    chamberShift           = cLoptions.GOoptions.chamberShift;
    
    LattStruct.LattData.geometry.ref.Magnets.HAperture  = magHAperture_ref*ones(npo_ref,1);
    LattStruct.LattData.geometry.ref.Chambers.HAperture = chamberHAperture_ref*ones(npo_ref,1);

    [x2d_magup_ref,y2d_magup_ref] = shiftPath(x2d_ref,y2d_ref,a2d_ref,...
                               magHAperture_ref*ones(npo_ref,1),1:npo_ref);   
    LattStruct.LattData.geometry.ref.Magnets.x2d_up=x2d_magup_ref;
    LattStruct.LattData.geometry.ref.Magnets.y2d_up=y2d_magup_ref;

    [x2d_magdown_ref,y2d_magdown_ref] = shiftPath(x2d_ref,y2d_ref,a2d_ref,...
                               -magHAperture_ref*ones(npo_ref,1),1:npo_ref);  
    LattStruct.LattData.geometry.ref.Magnets.x2d_down=x2d_magdown_ref;
    LattStruct.LattData.geometry.ref.Magnets.y2d_down=y2d_magdown_ref;
    
    [x2d_chaup_ref,y2d_chaup_ref] = shiftPath(x2d_ref,y2d_ref,a2d_ref,...
                               chamberHAperture_ref*ones(npo_ref,1),1:npo_ref);   
    LattStruct.LattData.geometry.ref.Chambers.x2d_up=x2d_chaup_ref;
    LattStruct.LattData.geometry.ref.Chambers.y2d_up=y2d_chaup_ref;

    [x2d_chadown_ref,y2d_chadown_ref] = shiftPath(x2d_ref,y2d_ref,a2d_ref,...
                               -chamberHAperture_ref*ones(npo_ref,1),1:npo_ref); 
    LattStruct.LattData.geometry.ref.Chambers.x2d_down=x2d_chadown_ref;
    LattStruct.LattData.geometry.ref.Chambers.y2d_down=y2d_chadown_ref;

    LattStruct.LattData.geometry.ref.chThick=distPath(x2d_chaup_ref,y2d_chaup_ref,x2d_chadown_ref,y2d_chadown_ref);
    switch cLoptions.GOoptions.GOmode
        case 1
            if (isfield(cLoptions,'RBfams'))
                if (not(isempty(cLoptions.RBfams)))
                    if (verboselevel>0)
                        fprintf('%s cLatt: calculating magnet and chamber geometry\n', datetime);
                    end
                    nRBs=numel(cLoptions.RBfams);
                    I_RBs=[];
                    for i=1:nRBs
                        Is = find(atgetcells(ACHRO_SP, 'FamName', cLoptions.RBfams{i}));
                        I_RBs=[I_RBs;Is];
                    end
                    dxrb = zeros(length(x2d),1);
                    for i=1:numel(I_RBs)
                        j=I_RBs(i);
                        RB=ACHRO_SP{j};
                        theta = RB.BendingAngle;
                        rho = RB.Length/theta;
                        dxrb(j) = 1/(RB.PolynomB(2)*rho);
                    end
                    [x2d_magce, y2d_magce] = shiftPath(x2d,y2d,a2d,...
                                             -dxrb,I_RBs);  
                    a2d_magce  = a2d;
                    dist_magceTobeam = distPath(x2d,y2d,x2d_magce,y2d_magce)./cos(a2d_magce');
                    magHAperture(1:npo) = magHAperture_ref;
                    magHAperture(I_RBs) = magHAperture_ref +...
                                          dist_magceTobeam(I_RBs) ; 


                    x2d_chace = x2d;
                    y2d_chace = y2d;
                    a2d_chace = a2d;
                    chaHAperture(1:npo) = chamberHAperture_ref;
                else
                    fprintf('%s cLatt Warning: Reverse Bend families not available for Magnet Centre calculation...\n', datetime);
                    x2d_magce  = x2d;
                    y2d_magce  = y2d;
                    dist_magceTobeam = zeros(npo,1);                                   
                    magHAperture(1:npo) = magHAperture_ref;
                    chaHAperture(1:npo) = chamberHAperture_ref;
                    x2d_chace = x2d;
                    y2d_chace = y2d;
                    a2d_chace = a2d;
                end
            else
                fprintf('%s cLatt Warning: Reverse Bend families not available for Magnet Centre calculation...\n', datetime);
                x2d_magce  = x2d;
                y2d_magce  = y2d;
                a2d_magce  = a2d;
                dist_magceTobeam = zeros(npo,1);                                   
                magHAperture(1:npo) = magHAperture_ref;
                chaHAperture(1:npo) = chamberHAperture_ref;
                x2d_chace = x2d;
                y2d_chace = y2d;
                a2d_chace = a2d;
            end

        case 2
            if (verboselevel>0)
                 fprintf('%s cLatt: calculating magnet and chamber geometry\n', datetime);
            end
            magHAperture(1:npo) = magHAperture_ref;
            x2d_magce = x2d;
            y2d_magce = y2d;
            a2d_magce = a2d;
            dist_magceTobeam = distPath(x2d,y2d,x2d_magce,y2d_magce)./cos(a2d_magce');                                   
            chaHAperture(1:npo) = chamberHAperture_ref;
            x2d_chace = x2d;
            y2d_chace = y2d;
            a2d_chace = a2d;
            
        case 3
            iBPM = findcells(ACHRO_ref_SP,'FamName','BPM');
            if (isempty(iBPM))
                iBPM=findcells(ACHRO_ref_SP,'FamName','mon');
            end
            %
            % The region of interest for the chamber shift is from the
            % third (U1/BPM-01) to the eighth (U5-BPM01) beam position
            % monitor. The chamber returns to nominal alomng the chamber
            % runnig from the second (M2/BPM-02) to the third BPMs and from
            % the eigth to the nineth (M2_BPM-01) BPMs
            %
            SBPM=findspos(ACHRO_ref_SP,iBPM);           
            chaHShift = zeros(npo_ref,1);
            chaHShift(iBPM(3):iBPM(8))=chamberShift;
 
            for i=iBPM(2):iBPM(3)
                chaHShift(i)=chamberShift/(SBPM(3)-SBPM(2))*(s2d_ref(i)-SBPM(2));
            end

            for i=iBPM(8):iBPM(9)
                chaHShift(i)=chamberShift-chamberShift/(SBPM(9)-SBPM(8))*(s2d_ref(i)-SBPM(8));
            end
 
            [x2d_chace,y2d_chace]=shiftPath(x2d_ref,y2d_ref,a2d_ref, chaHShift,1:npo_ref);
            a2d_chace = a2d_ref;
            chaHAperture(1:npo_ref) = chamberHAperture_ref;

            [x2d_magce,y2d_magce]=shiftPath(x2d_ref,y2d_ref,a2d_ref, chaHShift,1:npo_ref);
            
            y2d_magce = yint(x2d_magce,y2d_magce,x2d);
            a2d_magce = yint(x2d_magce,a2d_ref,x2d);
            x2d_magce = x2d;

            if (isfield(cLoptions,'RBfams'))
                if (not(isempty(cLoptions.RBfams)))

                    if (verboselevel>0)
                        fprintf('%s cLatt: calculating magnet and chamber geometry\n', datetime);
                    end
                    nRBs=numel(cLoptions.RBfams);
                    I_RBs=[];
                    for i=1:nRBs
                        Is = find(atgetcells(ACHRO_SP, 'FamName', cLoptions.RBfams{i}));
                        I_RBs=[I_RBs;Is];
                    end
                    dxrb = zeros(length(x2d),1);
                    for i=1:numel(I_RBs)
                        j=I_RBs(i);
                        RB=ACHRO_SP{j};
                        theta = RB.BendingAngle;
                        rho = RB.Length/theta;
                        dxrb(j) = 1/(RB.PolynomB(2)*rho) + chamberShift;
                    end
                    [x2d_magce, y2d_magce] = shiftPath(x2d_magce,y2d_magce,...
                                    a2d_magce,-dxrb,I_RBs);     
                    dist_magceTobeam   = distPath(x2d,y2d,x2d_magce,y2d_magce)./cos(a2d_magce');
                    magHAperture(1:npo) = magHAperture_ref;                
                    
                else
                    fprintf('%s cLatt Warning: Reverse Bend families not available for Magnet Centre calculation...\n', datetime);
                    x2d_magce  = x2d;
                    y2d_magce  = y2d;
                    dist_magce = zeros(length(x2d,1));
                    magH
                end
            else
                fprintf('%s cLatt Warning: Reverse Bend families not available for Magnet Centre calculation...\n', datetime);
                x2d_magce  = x2d;
                y2d_magce  = y2d;
                dist_magce = zeros(length(x2d,1));
            end            
        otherwise
                fprintf('%s cLatt Warning: Unknow GOmode %2i \n', datetime, cLoptions.GOmode);
                x2d_magce  = x2d;
                y2d_magce  = y2d;
                dist_magce = zeros(length(x2d,1));    
    end

    [x2d_magup,y2d_magup]     = shiftPath(x2d_magce,y2d_magce,a2d_magce,...
                                      magHAperture,1:numel(x2d_magce)); 

    [x2d_magdown,y2d_magdown] = shiftPath(x2d_magce,y2d_magce,a2d_magce,...
                                      -magHAperture,1:numel(x2d_magce));   

    [x2d_chaup,y2d_chaup]     = shiftPath(x2d_chace,y2d_chace,a2d_chace,...
                                      chaHAperture,1:numel(x2d_chace)); 

    [x2d_chadown,y2d_chadown] = shiftPath(x2d_chace,y2d_chace,a2d_chace,...
                                      -chaHAperture,1:numel(x2d_chace)); 


    dist_magup   = distPath(x2d_magup_ref,y2d_magup_ref,x2d_magup,y2d_magup);
    dist_magdown = distPath(x2d_magdown_ref,y2d_magdown_ref,x2d_magdown,y2d_magdown);

    dist_magce   = distPath(x2d_ref,y2d_ref,x2d_magce,y2d_magce);

    dist_chaup   = distPath(x2d_chaup_ref,y2d_chaup_ref,x2d_chaup,y2d_chaup);

    dist_chadown = distPath(x2d_chadown_ref,y2d_chadown_ref,x2d_chadown,y2d_chadown);
    
    dist_chace = distPath(x2d_ref,y2d_ref,x2d_chace,y2d_chace);
   
    dist_chaceTobeam = distPath(x2d,y2d,x2d_chace,y2d_chace)./cos(a2d_chace');

    

    dist_chaupTobeam   = distPath(x2d,y2d,x2d_chaup,y2d_chaup)./cos(a2d_chace');
    dist_chadownTobeam = distPath(x2d,y2d,x2d_chadown,y2d_chadown)./cos(a2d_chace');   
    effectiveAperture = min(dist_chaupTobeam ,dist_chadownTobeam);

    LattStruct.LattData.geometry.Magnets.Centre.x2d = x2d_magce;
    LattStruct.LattData.geometry.Magnets.Centre.y2d = y2d_magce;
    LattStruct.LattData.geometry.Magnets.Centre.Deviation = dist_magce;
    LattStruct.LattData.geometry.Magnets.HAperture = magHAperture;

    LattStruct.LattData.geometry.Magnets.Walls.x2d_up   = x2d_magup;
    LattStruct.LattData.geometry.Magnets.Walls.y2d_up   = y2d_magup;
    LattStruct.LattData.geometry.Magnets.Walls.x2d_down = x2d_magdown;
    LattStruct.LattData.geometry.Magnets.Walls.y2d_down = y2d_magdown;

    LattStruct.LattData.geometry.Magnets.Walls.dev_up   = dist_magup;
    LattStruct.LattData.geometry.Magnets.Walls.dev_down = dist_magdown;

    LattStruct.LattData.geometry.Magnets.magceTobeam = dist_magceTobeam;

    LattStruct.LattData.geometry.Chambers.Centre.x2d = x2d_chace;
    LattStruct.LattData.geometry.Chambers.Centre.y2d = y2d_chace;
    LattStruct.LattData.geometry.Chambers.Centre.Deviation = dist_chace;
    LattStruct.LattData.geometry.Chambers.HAperture = chaHAperture;

    LattStruct.LattData.geometry.Chambers.Walls.x2d_up   = x2d_chaup;
    LattStruct.LattData.geometry.Chambers.Walls.y2d_up   = y2d_chaup;
    LattStruct.LattData.geometry.Chambers.Walls.x2d_down = x2d_chadown;
    LattStruct.LattData.geometry.Chambers.Walls.y2d_down = y2d_chadown;

    LattStruct.LattData.geometry.Chambers.Walls.dev_up   = dist_chaup;
    LattStruct.LattData.geometry.Chambers.Walls.dev_down = dist_chadown;

    LattStruct.LattData.geometry.Chambers.chaceTobeam = dist_chaceTobeam;
    LattStruct.LattData.geometry.Chambers.effectiveAperture = effectiveAperture;
end
    
%% Calculates magnet challenge levels
if ((basicf||allf)&&not(isempty(fieldnames(MagnetStrengthLimits)))...
                  &&(not(isempty(cLoptions.eqfam)))...
                  &&(not(isempty(cLoptions.eqsca))))
    if (verboselevel>0)
        fprintf('%s cLatt: calculating magnet challenge levels \n', datetime);
    end
    eqfam = cLoptions.eqfam;
    eqsca = cLoptions.eqsca;
    All_famsO = cLoptions.All_famsO;
    CLv = chalevel(MagnetStrengthLimits,'mode','LS',...
                   'ACHRO',ACHRO,'eqfam', eqfam,'eqsca',...
                    eqsca,'Fams',All_famsO,'verbose',verboselevel-1)';  
    LattStruct.LattData.CLv=CLv;
    if (isfield(CLv,'outputs'))
        XAllO=CLv.outputs.X0;
        LattStruct.LattData.XAllO=XAllO;
    end
end

%% Creates full ring structure if not yet available
if (isempty(RING))
     if (verboselevel>0)
        fprintf('%s cLatt Warning: Full ring input not available, creating ring from achromat...\n',datetime);
     end
     RING = atenable_6d(SetBPMWeights(achromat2ring(ACHRO)));
end
cavpts  = find(atgetcells(RING, 'FamName', 'CAVITY'));
if (isempty(cavpts))
    cavpts  = find(atgetcells(RING, 'FamName', 'CAV'));
end
if (isempty(cavpts))
    cavpts  = find(atgetcells(RING, 'FamName', 'cav'));
end
if (not(isempty(cavpts)))
   RING = atSetRingProperties(RING,'rf_frequency','nominal',...
          'cavpts',cavpts,'rf_voltage',V0);
   RINGGRD=calculateGirderMaps(RING);
   LattStruct.LattData.RINGGRD=RINGGRD;
else
   fprintf('%s cLatt Warning: Full ring input does not contain cavities...\n',datetime);
end

if(isempty(RINGGRD))
    if (verboselevel>0)
        fprintf('%s cLatt Warning: full ring cell array not available\n', datetime);
    end
end

%% Calculates field and gradient profiles
if (basicf||allf)
    if (verboselevel>0)
        fprintf('%s cLatt: calculating field and gradient profiles \n', datetime);
    end
    FG = calcFields(ACHRO,cLoptions.All_famsO,'desc',lattname,'split',split);
    LattStruct.LattData.FG=FG;
end

%% Calculates atsummary for full ring (or an achromat if ring not available)
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
%% Generates lattice at zero chromaticity
if ( (basicf||allf)&&not(isempty(cLoptions.chrom_fams)))
    if (verboselevel>0)
         fprintf('%s cLatt: creating zero chromaticity lattice\n',datetime);
    end
    TolChrom     = cLoptions.DAoptions.TolChrom; % Chromaticity tolerances
    Nitchro      = cLoptions.DAoptions.Nitchro;  % Max n. iterations of chromaticity correction
    chrom_fams   = cLoptions.chrom_fams;
    try 
        [ACHRO_zc, ~, ~]=fitchroit(ACHRO, chrom_fams, [0 0], Nitchro, TolChrom); 
    catch ME
        fprintf('%s cLatt: Error in fitting zero chromaticity: ', datetime);
        fprintf('Error message was:%s \n',ME.message);
        ACHRO_zc={};
    end
    LattStruct.LattData.ACHROMAT_ZC = ACHRO_zc;
end
ACHRO_zc=LattStruct.LattData.ACHROMAT_ZC;
%% Generates lattice with all sextupoles at zero
if ( (basicf||allf)&&not(isempty(cLoptions.sext_fams)))
    if (verboselevel>0)
         fprintf('%s cLatt: creating lattice with all sextupoles at zero\n',datetime);
    end
    SFams=cLoptions.sext_fams;
    nSFams = numel(SFams);
    ACHRO_zs = ACHRO;
    for i=1:nSFams
        I_sext   = find(atgetcells(ACHRO, 'FamName', SFams{i}));
        ACHRO_zs = atsetfieldvalues(ACHRO_zs, I_sext, 'PolynomB',{1,3}, 0.0); 
    end
    LattStruct.LattData.ACHROMAT_zs = ACHRO_zs;
    ats_zs=atsummary(ACHRO_zs);
    LattStruct.LattPerf.atsummary.achrnatchrom=ats_zs.chromaticity;
    LattStruct.LattPerf.atsummary.achrtunes=ats_zs.Itunes;
end

% waitbar
frac=frac+dfrac;
waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
if getappdata(fb,'canceling')
       fprintf ('Cancelling at %5.2f %% \n', frac);
       delete(fb);
       exitflag = 'cancelled';
       return
end
%% Evaluates DAs for full RING
if (not(isempty(RINGGRD)))
%% Calculates Dynamic aperture for full ring without errors on (x,y) plane
  if (DAxyf||allf||DAsf)
    if (verboselevel>0)
      fprintf('%s DA calculation without errors: on-momentum. \n', datetime);
    end
    DAS_0 = calcDA(RINGGRD,cLoptions.DAoptions, 'mode', 'xy', 'dp', 0.00, 'verbose', verboselevel-1);
    LattStruct.LattPerf.DA.xy_0=DAS_0;
% waitbar
    frac=frac+dfrac;
    waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
    if getappdata(fb,'canceling')
       fprintf ('Cancelling at %5.2f %% \n', frac);
       delete(fb);
       exitflag = 'cancelled';
       return
    end
    if (verboselevel>0)
      fprintf('%s DA calculation without errors: +3 %%. \n', datetime);
    end
    DAS_p3 = calcDA(RINGGRD,cLoptions.DAoptions,'mode', 'xy', 'dp',+0.03, 'verbose', verboselevel-1);
    LattStruct.LattPerf.DA.xy_p3=DAS_p3; 
% waitbar
    frac=frac+dfrac;
    waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
    if getappdata(fb,'canceling')
       fprintf ('Cancelling at %5.2f %% \n', frac);
       delete(fb);
       exitflag = 'cancelled';
       return
    end
    if (verboselevel>0)
      fprintf('%s DA calculation without errors: -3 %%. \n', datetime);
    end
    DAS_m3 = calcDA(RINGGRD,cLoptions.DAoptions,'mode', 'xy', 'dp',-0.03, 'verbose', verboselevel-1);
    LattStruct.LattPerf.DA.xy_m3=DAS_m3;
% waitbar
    frac=frac+dfrac;
    waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
    if getappdata(fb,'canceling')
       fprintf ('Cancelling at %5.2f %% \n', frac);
       delete(fb);
       exitflag = 'cancelled';
       return
    end
  else
      frac=frac+3*dfrac;
  end
  
%% Calculates Dynamic aperture for full ring without errors on (x,dp) and (y,dp) planes
  if (DAxydpf||allf||DAsf)
    if (verboselevel>0)
      fprintf('%s DA calculation without errors: (xy,dp) planes \n', datetime);
    end
    DAS = calcDA(RINGGRD,cLoptions.DAoptions, 'mode', 'xydp','verbose', verboselevel-1);
    LattStruct.LattPerf.DA.xydp=DAS;
  end
% waitbar
  frac=frac+dfrac;
  waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
  if getappdata(fb,'canceling')
      fprintf ('Cancelling at %5.2f %% \n', frac);
      delete(fb);
      exitflag = 'cancelled';
      return
  end  
%% Calculates Dynamic aperture for full ring with errors on xy plane
  if ((DAdistxyf||allf||DAsf)&&(not(isempty(cLoptions.ringtune_fams))))
    if (verboselevel>0)
      fprintf('%s DA calculation with errors: on-momentum. \n', datetime);
    end
    DAdist = calcDAdist(RINGGRD,cLoptions.ErrorModel,cLoptions.DAoptions,...
             'tunfams',cLoptions.ringtune_fams,'mode','xy',...
             'corrorb',cLoptions.corrorb,'corrtun',cLoptions.corrtun,...
             'frac',cLoptions.tunfrac,'nseeds',cLoptions.nseeds,...
             'verbose', verboselevel-1);
    LattStruct.LattPerf.DAdist.xy = DAdist;
  end
  
% waitbar
  frac=frac+dfrac;
  waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
  if getappdata(fb,'canceling')
      fprintf ('Cancelling at %5.2f %% \n', frac);
      delete(fb);
      exitflag = 'cancelled';
      return
  end

%% Calculates Dynamic aperture for full ring with errors on xdp an ydp planes
  if ((DAdistxydpf||allf||DAsf)&&(not(isempty(cLoptions.ringtune_fams))))
    if (verboselevel>0)
      fprintf('%s DA calculation with errors: off-momentum. \n', datetime);
    end
    DAdist = calcDAdist(RINGGRD,cLoptions.ErrorModel,cLoptions.DAoptions,...
             'tunfams',cLoptions.ringtune_fams,'mode','xydp',...
             'corrorb',cLoptions.corrorb,'corrtun',cLoptions.corrtun,...
             'frac',cLoptions.tunfrac,'nseeds',cLoptions.nseeds,...
             'verbose', verboselevel-1);
    LattStruct.LattPerf.DAdist.xydp = DAdist;
    % waitbar
    frac=frac+dfrac;
    waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
    if getappdata(fb,'canceling')
      fprintf ('Cancelling at %5.2f %% \n', frac);
      delete(fb);
      exitflag = 'cancelled';
      return
    end
  end
else
    if (verboselevel>0)
        fprintf('%s cLatt Warning : RINGRD structure not available. \n', datetime);
    end
end

%% Calculates tune maps for an achromat
if (not(isempty(ACHRO_zc)))
  TMoptions=cLoptions.TMoptions;  
%% Tune map along x and y axes
  if(TM_xyf||allf||TMsf)
      if (verboselevel>0)
        fprintf('%s Tune Map xy (ADTS) \n', datetime);
      end
      
      TMoptions.mode='xy';
      tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
                'plottype','xy');
      LattStruct.LattPerf.TM.xy=tunemap;
  end

%% Tune map on a grid of points in (x,y) plane.
  if (TM_gridxyf||allf||TMsf)
       if (verboselevel>0)
            fprintf('%s Tune Map grid (x,y) (ADTS). \n', datetime);
       end
       
       TMoptions.mode='gridxy';
       tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
                 'plottype','gridxy');
       LattStruct.LattPerf.TM.gridxy=tunemap;
  end
%% Tune map on a grid of points in (x,dp) plane
  if (TM_gridxdpf||allf||TMsf)
    if (verboselevel>0)
        fprintf('%s Tune Map grid (x,dp). \n', datetime);
    end

    TMoptions.mode='gridxdp';
    tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
               'plottype','gridxdp');
    LattStruct.LattPerf.TM.gridxdp=tunemap;
  end

%% Tune map on a grid of points in (y,dp) plane
  if (TM_gridydpf||allf||TMsf)
    if (verboselevel>0)
        fprintf('%s Tune Map grid (y,dp). \n', datetime);
    end
    
    TMoptions.mode='gridydp';
    tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
               'plottype','gridydp');
    LattStruct.LattPerf.TM.gridydp=tunemap;
  end

%% Tune diffusion map on a grid of points in (x,y) plane
  if (TM_difxyf||allf||TMsf)
    if (verboselevel>0)
        fprintf('%s Tune diffusion map (x,y). \n', datetime);
    end

    TMoptions.mode='difxy';
    tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
               'plottype','difxy');
    LattStruct.LattPerf.TM.difxy=tunemap;
  end

%% Tune diffusion map on a grid of points in (x,dp) plane
  if (TM_difxdpf||allf||TMsf)
    if (verboselevel>0)
        fprintf('%s Tune diffusion map (x,dp). \n', datetime);
    end

    TMoptions.mode='difxdp';
    tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
               'plottype','difxdp');
    LattStruct.LattPerf.TM.difxdp=tunemap;
  end

%% Tune diffusion map on a grid of points in (y,dp) plane
  if (TM_difydpf||allf||TMsf)
    if (verboselevel>0)
        fprintf('%s Tune diffusion map (y,dp). \n', datetime);
    end

    TMoptions.mode='difydp';
    tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
               'plottype','difydp');
    LattStruct.LattPerf.TM.difydp=tunemap;
  end

%% Chromatic tune map
  if (TM_chrof||allf||TMsf)
    if (verboselevel>0)
        fprintf('%s Chromatic tune map. \n', datetime);
    end

    TMoptions.mode='chro';
    tunemap = calcTuneMap(ACHRO_zc,TMoptions,...
               'plottype','chro');
    LattStruct.LattPerf.TM.chro=tunemap;
  end
else
   fprintf('%s cLatt Warning: ACHRO_zc (zero chromaticty achromat) structure not available. \n', datetime);
end
% waitbar
frac=frac+dfrac;
waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
if getappdata(fb,'canceling')
   fprintf ('Cancelling at %5.2f %% \n', frac);
   delete(fb);
   exitflag = 'cancelled';
   return
end
%% Calculates tune map for a full ring with errors
if (not(isempty(RINGGRD)))
    if (TM_distf||allf)
        TMdist=calcTMdist(RINGGRD,cLoptions.ErrorModel,...
            cLoptions.TMoptions,'tunfams',cLoptions.ringtune_fams,...
            'nseeds',1,'mode','difxy','verbose',verboselevel-1);
        LattStruct.LattPerf.TMdist=TMdist;
    end
end
% waitbar
frac=frac+dfrac;
waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
if getappdata(fb,'canceling')
   fprintf ('Cancelling at %5.2f %% \n', frac);
   delete(fb);
   exitflag = 'cancelled';
   return
end

%% Calculates LMA for full ring without errors

if (not(isempty(RINGGRD)))
    if (LMAf||allf)
        if (verboselevel>0)
            fprintf('%s Local Momentum Aperture without errors \n', datetime);
        end
        LMA=calcLMA(RINGGRD,cLoptions.MAoptions,'verbose',verboselevel-1);
        LattStruct.LattPerf.LMA=LMA;
    end
else
    if (verboselevel>0)
        fprintf('%s cLatt Warning: RING structure not available for LMA calculation. \n', datetime);
    end
end
% waitbar
  frac=frac+dfrac;
  waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
  if getappdata(fb,'canceling')
      fprintf ('Cancelling at %5.2f %% \n', frac);
      delete(fb);
      exitflag = 'cancelled';
      return
  end
%% Calculates LMA for a full ring with errors
if (not(isempty(RINGGRD)))
    if ((LMAdistf||allf)&&(not(isempty(cLoptions.ringtune_fams))))
        if (verboselevel>0)
            fprintf('%s Local Momentum Aperture with errors \n', datetime);
        end
        LMAdist=calcLMAdist(RINGGRD,cLoptions.ErrorModel,...
             cLoptions.MAoptions,'verbose',verboselevel-1,...
             'corrorb',cLoptions.corrorb,'corrtun',cLoptions.corrtun,...
             'tunfams',cLoptions.ringtune_fams,'nseeds',cLoptions.nseeds);
        LattStruct.LattPerf.LMAdist=LMAdist;
    end
else
    if (verboselevel>0)
        fprintf('%s Error: RING structure not available for LMAdist calculation. \n', datetime);
    end
end
% waitbar
  frac=frac+dfrac;
  waitbar(frac/100,fb,strcat(sprintf('%3.0f %s',frac,'%')));
  if getappdata(fb,'canceling')
      fprintf ('Cancelling at %5.2f %% \n', frac);
      delete(fb);
      exitflag = 'cancelled';
      return
  end
%% Calculates Touschek lifetime for a ring without errors
if (not(isempty(RINGGRD)))
    if (TLTf||allf)
        if (verboselevel>0)
            fprintf('%s Touschek lifetime without errors \n', datetime);
        end
        TL=calcTLT(RINGGRD,cLoptions.TLoptions,cLoptions.MAoptions,...
            'LMA',LattStruct.LattPerf.LMA,...
            'LMAPeriods',1,'verbose',verboselevel-1);
        LattStruct.LattPerf.TL=TL;
    end
else
    if (verboselevel>0)
        fprintf('%s cLatt Warning: RING structure not available for Touschek lifetime calculation. \n', datetime);
    end
end
%% Calculates Touschek lifetime for a full ring with errors
%
if (not(isempty(RINGGRD)))
    if ((TLTdistf||allf)&&(not(isempty(cLoptions.ringtune_fams))))
         if (verboselevel>0)
            fprintf('%s Touchek lifetime with errors \n', datetime);
        end
        TLdist=calcTLTdist(RINGGRD,cLoptions.ErrorModel,...
           cLoptions.TLoptions,cLoptions.MAoptions,...
           'LMAdist',LattStruct.LattPerf.LMAdist,...
           'verbose',verboselevel-1,'corrorb',cLoptions.corrorb,...
           'corrtun',cLoptions.corrtun,'tunfams',cLoptions.ringtune_fams,...
           'nseeds',cLoptions.nseeds);
           LattStruct.LattPerf.TLdist=TLdist;
    end
else
    if (verboselevel>0)
        fprintf('%s Error: RING structure not available for TLTdist calculation. \n', datetime);
    end
end

%% Clean up and exit
delete(fb);
exitflag = 'normal';
fprintf('%s Lattice structure creation/update/evaluation completed. \n', datetime);
fprintf(' ************* \n');

%% Auxiliary functions

function ringW=SetBPMWeights(ring)
% In order to avoid issues with orbit correction routines, which all seem
% to assume that we have as many vertical correctors as we do BPMs, we need
% to add weight information to the BPMs (and make sure correction routines
% make use of it).
%
    ringW=ring;
    I = findcells(ringW,'FamName','BPM');
    if (isempty(I))
        I=findcells(ringW,'FamName','mon');
    end
    
    if (not(isempty(I)))
        ringW{I(1)}.Weight = [0 0]; I = I(2:end);
        for n = 1:numel(I)
            if mod(n-2,10) == 0
                ringW{I(n)}.Weight = [1 1e-3];
            else
                ringW{I(n)}.Weight = [1 1];
            end
        end
    else
        fprintf('%s SetBPMWeigths: Warning - no BPMs found. \n', datetime);
    end

function [x2d_s,y2d_s]=shiftPath(x2d,y2d,a2d,shift,Is)
%
% returns a path parallel to the input path shifted by a 
% value dependent on the point along the path
    x2d_s=x2d;
    y2d_s=y2d;

    for j = 1:numel(Is)
        i=Is(j);
        x2d_s(i)   = x2d(i) - shift(i)*sin(a2d(i));
        y2d_s(i)   = y2d(i) + shift(i)*cos(a2d(i));
    end

function dist = distPath(x2d,y2d,x2d_s,y2d_s)
%
% returns vertical distance between two paths
%
    dist = zeros(length(x2d_s),1);
    [~, ia, ~] = unique(x2d);
    x2d_u = x2d(ia);
    y2d_u = y2d(ia);
    for i = 1:length(x2d_s)
        yinter_s=interp1(x2d_u,y2d_u,x2d_s(i));
        dist(i)    = abs(yinter_s - y2d_s(i));
    end

function yint = yint(x2d,y2d,x2dint)
%
    yint=zeros(length(x2dint),1);
    [~, ia, ~] = unique(x2d);
    x2d_u = x2d(ia);
    y2d_u = y2d(ia);
    for i=1:length(x2dint)
        yint(i)=interp1(x2d_u,y2d_u,x2dint(i));
    end

              