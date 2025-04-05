%% m4Uc_Template
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
% 2024/07/17 : incrporrated Saroj's mod_IntSteps function
% 2024/07/19 : incorporated geometry analysis and removed diary logs
% 2024/08/07: added the input of bucket height and harmonic number 

%% Lattice specific data for ib03-03-09-08
%% leading idea
%% start from 241106_b03_09_08 and gradually move to i-type lattice,
%% preserving phase advances between sextupoles 
%% first step, reproducing the old result with the new release of MAXIV_addition

load 'm4U_250404_ib03_03_09_08.mat';
lattname= 'm4Uc-250404-ib03-03-09-08'; 

desc = 'stems from 241106_b03_09_08 --> moving to i-type case';

ACHRO = ACHROMAT; % 

%load CandidateLattices/m4_standard/m4_standard.mat
% load ../m4_standard/m4_standard.mat
% load MagnetStrengthLimits.mat

load /gpfs/offline1/staff/common/marapo/MAXIV/DESIGN/Lattices/MAXIV/NewLattices/SCRIPTS/MAX4U_PFT/at-AT2_PFT-MAX4U/MAX4U/CandidateLattices/m4_standard/m4_standard.mat
load /gpfs/offline1/staff/common/marapo/MAXIV/DESIGN/Lattices/MAXIV/NewLattices/SCRIPTS/MAX4U_PFT/at-AT2_PFT-MAX4U/MAX4U/MagnetStrengthLimits.mat



% Initialize physical apertures
for i=1:length(ACHRO)
    ACHRO{i}.EApertures=[0.011 0.011];

end

clear('cLoptions');

cLoptions.All_famsO={}; %   optional, If empty m4_cLatt finds out the magnet 
%                           family names. Use this in case a specific order
%                           of family names is desired.

cLoptions.ringtune_fams = {'Q1';'Q2'};   % magnet families for tune matching
cLoptions.chrom_fams    = {'S1','S5'}; % magnet families for chromaticity matching
cLoptions.DAoptions.chroms0 = [2 2]/20;  
cLoptions.sext_fams     = {'S1';'S2';'S3';'S4';'S5'}; % list of all sextupole families 
cLoptions.RBfams        = {'R1'}; % reverse bend families
cLoptions.corchro       = false;% false;
cLoptions.corrorbf      = true;
cLoptions.corrtunf      = true;

cLoptions.DAoptions.DAmode        = 'smart_out'; 
cLoptions.DAoptions.nturns        = nan;
cLoptions.DAoptions.smart_divider = 1.4; 
cLoptions.DAoptions.da_yclip      = 1e+19;

cLoptions.OCoptions.cflags        = [true false];

cLoptions.ErrorModel = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);

cLoptions.eqfam = {'dipm';'dip';'dip';...                  
                   'Qf_Qfm';'Oxx_Oxy';'Oxx_Oxy';'Oyy';...
                   'Qfend_Qdend';'Qfend_Qdend';'Qf_Qfm';...
                   'Qf_Qfm';'Sdend';'Sfm';'Sfi_Sfo';'Sfi_Sfo';'Sfi_Sfo'};
cLoptions.eqsca = [1 1 (3.0+0.88)/3.0 1  1 1 1 1 1 1 1 1 1 1 1 1];

cLoptions.sumTabfams = {'D1';'D2';'D3';...
                        'Q1';'Q2';'Q3';'Q4';'';'';'';...
                        'R1';'';'';...
                        'S1';'S2';'S3';'S4';'S5';'';...
                        'O1';'O2';'O3';...
                        '';'';...
                        '';'';''...
                        };

cLoptions.sumTabsca  = [2 4 1 ...
                        1 1 1 1 0 0 0 ...
                        1 0 0 ...
                        1 1 1 1 1 0 ...
                        1 1 2 ...
                        0 0 ...
                        0 0 0];

cLoptions.eqsumTabfams = {'dipm';'dip';'dip';...
                          'Qfend';'Qdend';'Qfm_1';'Qf_U2_1';'';'';'';...
                          'Rb_U3_1';'';'';...
                          '';'Sfm';'';'Sfo';'Sfi';'';...
                          'Oxx';'Oxy';'Oyy';...
                          '';''; ...
                          'Sdendcomb';'Sdcomb';'Sdcomb'...
                          };


% V0 = 1.8E6;
bh=0.1;   % Bucket Height
harm=176; % Harmonic number 
corchrof=0; 
%% Run cLatt options
V0 = 1.0e6; 
m4UT = m4Uc_Latt(ACHRO,lattname,desc,cLoptions,m4_standard.ACHROMAT,MagnetStrengthLimits,'bh',bh,'V0',V0,'harm',harm,'corchro',corchrof);


% Below an example of how the 'plotLatt' function can be used to prodcue
% plots from the m4UT structure and save the results on a file.
% plotLatt(m4UT,'all','ymaxplot_dm',0.004,'zoom',2.0,'ymaxplot',0.004,'xminplot',-0.010,'xmaxplot',0.01,'dpminplotLMA',-0.25,'dpmaxplotLMA',0.25,'nogrid','save');