%% m4Uc_240701_f01_01_01_01
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
% A quick look at the m4Uc_Latt function reveals the syntx of calls to
%  the cLatt function
%

%% History
% 2024/07/07 : first verstion
% 2024/07/10 : break up of the calculations into two steps - with and
%              without errors
% 2024/07/13 : restructured to call m4Uc_Latt as a function
%% Lattice specific data
lattname = 'm4U-240807-f01-03-01-01';
desc = 'Downloaded from JB tema folder on 20240807';

ACHRO = max_4u_f_0_20240807(); % The cell array with the AT2 lattice to be evaluated

% Initialize physical apertures
 for i=1:length(ACHRO)
    ACHRO{i}.EApertures=[0.011 0.011];
 end

clear('cLoptions');
cLoptions.All_famsO={}; %   optional, If empty m4_cLatt finds out the magnet 
%                           family names. Use this in case a specific order
%                           of family names is desired.
cLoptions.ringtune_fams = {'q1';'q2'};   % magnet families for tune matching
cLoptions.chrom_fams    = {'s2';'s3'}; % magnet families for chromaticity matching
cLoptions.sext_fams     = {'s1';'s2';'s3';'s4'}; % list of all sextupole families 
cLoptions.eqfam = {'dipm';'dipm';'dipm';'dipm';'dipm';'dipm';'dipm';...
                   'dipm';'dipm';'dipm';'dipm';'dipm';'dip';'dip';...
                   'dip';'dip';'dip';'dip';...
                   'Qf_Qfm';'Oxx_Oxy'; 'Oxx_Oxy'; 'Oyy';...
                   'Qfend_Qdend';'Qfend_Qdend';'Qf_Qfm';...
                   'Sdend'; 'Sfi_Sfo'; 'Sfi_Sfo';'Sfi_Sfo'};
cLoptions.RBfam = {'r1'};

cLoptions.eqsca = ones(29,1);


load(strcat(erase(atroot,'atmat'),'/MAX4U/MagnetStrengthLimits.mat'));
load(strcat(erase(atroot,'atmat'),'/MAX4U/CandidateLattices/m4_standard/m4_standard.mat'));
ACHRO_ref = m4_standard.ACHROMAT;
V0=1.8E6;
bh='auto';

%% Run cLatt options

m4UT = m4Uc_Latt(ACHRO,lattname,desc,cLoptions,ACHRO_ref,MagnetStrengthLimits,'V0',V0,'bh',bh,'corchro',true,'basonly');

%plotLatt(m4UT,'all','ymaxplot_dm',0.007,'zoom',2.0,'ymaxplot',0.007,'xminplot',-0.012,'xmaxplot',0.012,'dpminplotLMA',-0.20,'dpmaxplotLMA',0.20,'nogrid','xmaxplot_dm',0.012,'xminplot_dm',-0.012,'caxrange',[-10 0],'caxrange_r',[-10 -5]);