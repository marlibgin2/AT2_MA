%% m4Uc_240316_b03_01_03_01
% This script is a template to set up workspace variables needed
% to run the "m4U_cLatt" script, which in turn, runs all "cLatt" function 
% options in a series of steps creating/overwriting a structure named
% "m4UT" in the matlab workspace, saving this variable
% in a "m4UT.mat" file and saving all matlab output in a logfile
%
% The intended workflow is 
%  1. Edit the statements in the "Lattice Specific Data" section below.
%  2. Make sure the variables ACHROGRD_a1 and MagnetStrengthLimits are 
%     available in the workspace (see below for their content)
%  3. Run the script.
%  4. Rename or copy the resulting "m4UT" structure accordingly.
%  5. Execution may be interrupted at any moment and the most recent 
%     m4UT scruture can be recovered by loading it from the file 
%     saved on disk. From that point on, further calculation and 
%     updates to the m4UT structure may be carried out by running 
%     cLatt with specific input options. See the m4Uc_Latt script for 
%     many possible calls to cLatt.
%  
% The structure "cLoptions" is stored as a field in the "m4UT" structure 
% and contains various optionallattice evaluation settings. These are
% listed in the "General initialisation" section of the m4Uc_Latt function 
% (which sets their default values). If values different from the defauls 
% are desired for fields in the cLoptions structure, these can be changed 
% at anytime, e.g., by directly editing the correspondig fields or by 
% typing:
%    >>cLpptions = m4UT.cLoptions; % extracts cLoptions to the workspace.
%    >>cLoptions.DAoptions.nturns = 2028 % implementi the desired changes -
%    in this example, we change the number of turns used in dynamic
%    aperture calcualations
%
%    >>m4UT=cLatt(m4UT, 'cLoptions',cLoptions,'verbose',1); % stores the
%                                          changes in the "m4UT" structure.
%    Note that such changes do not automatically re-run any evaluations 
%    that are affected by them - these need to be rerun explicitly 
%    to update the corresponding fields in the m4UT structure, e.g.
%
%    >> m4UT = cLatt(m4UT,'DAxy','verbose',2);
%    
% A quick look at the m4Uc_Latt function reveals the syntax of calls to
%  the cLatt function

%% History
% 2024/07/07 : first verstion
% 2024/07/10 : break up of the calculations into two setps - with and
%              without errors
% 2024/07/13 : restructured to call m4Uc_Latt as a function
%
%% Lattice specific data
lattname = 'm4U-240316-b03-01-03-01';
desc = 'MOGA_20240313T091640, Ind 12, Disp Matched, 3.49mrad RB';

m4U('b3','Full','rbk',3.49E-3);
DAoptions = LatticeOptData.DAoptions;
DAoptions.xmaxdas  =  0.010;
DAoptions.xmindas  = -0.010;
DAoptions.ymnaxdas =  0.007;


load('/home/pedtav/Documents/Max4U/MOGA_Scans/MOGA_20240313T091640.mat');


rp=ExMOGA(MOGAResults,12,'verbose','fitdisp','fitchrom','fittune',...
          'tunes',[55.15 16.20],'LatticeOptData',LatticeOptData,...
          'DAoptions',DAoptions);


ACHRO = rp.outputs.ACHROGRD; % The cell array with the AT2 lattice to be evaluated

% Initialize physical apertures
%for i=1:length(ACHRO)
%    ACHRO{i}.EApertures=[0.011 0.011];
%end

clear('cLoptions');
cLoptions.All_famsO=LatticeOptData.All_famsO; %   optional, If empty m4_cLatt finds out the magnet 
%                           family names. Use this in case a specific order
%                           of family names is desired.
cLoptions.ringtune_fams = LatticeOptData.ringtune_fams;   % magnet families for tune matching
cLoptions.chrom_fams    = LatticeOptData.chrom_fams; % magnet families for chromaticity matching
cLoptions.RBfams = LatticeOptData.RBfams;
cLoptions.eqfam = LatticeOptData.eqfam;
cLoptions.eqsca = LatticeOptData.eqsca;
cLoptions.ErrorModel = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);

load('/home/pedtav/Documents/Codes/AT/AT2.0/MAX4U/MagnetStrengthLimits.mat');
load('/home/pedtav/Documents/Codes/AT/AT2.0/MAX4U/CandidateLattices/m4_standard/m4_standard.mat');
ACHRO_ref = m4_standard.ACHROMAT;

%% Run cLatt options
m4UT = m4Uc_Latt(ACHRO,lattname,desc,cLoptions,ACHRO_ref,MagnetStrengthLimits);

%plotLatt(m4UT,'all','ymaxplot_dm',0.004,'zoom',2.0,'ymaxplot',0.004,'xminplot',-0.010,'xmaxplot',0.01,'dpminplotLMA',-0.25,'dpmaxplotLMA',0.25,'nogrid','save');
