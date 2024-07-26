%% m4Uc_jb5
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
% 2024/07/10 : break up of the calculations into two setps - with and
%              without errors
% 2024/07/13 : restructured to call m4Uc_Latt as a function
% 2024/07/19 : incorporated geometry analysis and removed diary logs
%% Lattice specific data
lattname = 'm4U-jb5';
desc = 'From Johan 20240718';

ACHRO = max_4u_sp_jb_5(); % The cell array with the AT2 lattice to be evaluated

% ACHRO = mod_IntSteps(ACHRO,10,10,3,1); % If necessary, fix missing NumIntSteps

% Initialize physical apertures
for i=1:length(ACHRO)
   ACHRO{i}.EApertures=[0.011 0.011];
end

clear('cLoptions');
cLoptions.All_famsO={}; %   optional, If empty m4_cLatt finds out the magnet 
%                           family names. Use this in case a specific order
%                           of family names is desired.
cLoptions.ringtune_fams = {'qf2';'qd'};   % magnet families for tune matching
cLoptions.chrom_fams    = {'sd1';'sf_h'}; % magnet families for chromaticity matching
cLoptions.sext_fams     = {'sd1';'sf_h'}; % list of all sextupole families 
cLoptions.eqfam = {'dip';'';'';'';'';'';'dip';...
                   'dip';'dip';'dip';'dip';'dip';'dip';'dip';...
                   'dip';'dip';'dip';'dip';...
                   'Qfend_Qdend';'Qf_Qfm';'Qf_Qfm';...
                   'Qfend_Qdend'; 'Oxx_Oxy'; 'Oxx_Oxy'; 'Oyy'; 'Sd';...
                   'Sdend'; 'Sfi_Sfo'; 'Sfi_Sfo'};
cLoptions.RBfams = {'qf1';'qf1_e'};

cLoptions.eqsca = ones(29,1);

cLoptions.GOoptions.GOmode = 1;
cLoptions.GOoptions.chamberHAperture   = 11.0E-3;
cLoptions.GOoptions.chamberTomagnetGap =  0.5E-3;
cLoptions.GOoptions.chamberThickness   =  1.0E-3;
cLoptions.GOoptions.chamberShift       =  4.0E-3;

cLoptions.ErrorModel = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);

load(strcat(erase(atroot,'atmat'),'/MAX4U/MagnetStrengthLimits.mat'));
load(strcat(erase(atroot,'atmat'),'/MAX4U/CandidateLattices/m4_standard/m4_standard.mat'));
ACHRO_ref = m4_standard.ACHROMAT;

%% Run cLatt options
m4UT = m4Uc_Latt(ACHRO,lattname,desc,cLoptions,ACHROGRD_a1,MagnetStrengthLimits);

%plotLatt(m4UT,'all','ymaxplot_dm',0.005,'zoom',2.0,'ymaxplot',0.005,'xminplot',-0.006,'xmaxplot',0.006,'dpminplotLMA',-0.10,'dpmaxplotLMA',0.10,'nogrid','xmaxplot_dm',0.006,'xminplot_dm',-0.006,'caxrange',[-10 0],'caxrange_r',[-10 -5],'save');