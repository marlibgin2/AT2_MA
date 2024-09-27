%% m4Uc_Template
% This script is a template to set up variables needed
% to run the "m4U_cLatt" function, which in turn, runs all "cLatt" function 
% options in a series of steps creating/overwriting a structure named
% "m4UT" in the matlab workspace, saving this variable
% in a "m4UT.mat" file and saving all matlab output in a logfile
%
% The intended workflow is 
%  1. Edit the statements in the "Lattice Specific Data" section below.
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
% and contains various optional lattice evaluation settings. These are
% listed in the "General initialisation" section of the m4Uc_Latt function 
% (which sets their default values). If values different from the defaults 
% are desired for fields in the cLoptions structure, these can be changed 
% at anytime, e.g., by directly editing the correspondig fields or by 
% typing:
%    >>cLpptions = m4UT.cLoptions; % extracts cLoptions to the workspace.
%    >>cLoptions.DAoptions.nturns = 2028 % implementi the desired changes -
%    in this example, we change the number of turns used in dynamic
%    aperture calculations
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
% 2024/07/17 : incorporated Saroj's mod_IntSteps function
% 2024/07/19 : incorporated geometry analysis and removed diary logs
%
%% Lattice specific data
lattname = 'm4U-240618-f02-02-01-01';
desc = 'From Åke, symmetrized and translated by Johan B. to Tracy 2024/07/02';

ACHRO = m4_20240618_M1a_QFAQFBQFCp175c_tracy(); % run or load tha ACHROMAT 
                                                % cell array containing the 
                                                % AT2 lattice for a
                                                % single achromat
                                                

% ACHRO = mod_IntSteps(ACHRO,10,10,3,1); % If necessary, fix missing NumIntSteps

% Initialize physical apertures
for i=1:length(ACHRO)
    ACHRO{i}.EApertures=[0.011 0.011];
end

clear('cLoptions');
cLoptions.All_famsO={}; %   optional, If empty, m4_cLatt finds out the magnet 
%                           family names. Use this in case a specific order
%                           of family names is desired.

cLoptions.ringtune_fams = {'qfend';'qdend'};   % magnet families for tune matching
cLoptions.chrom_fams    = {'sdqd','sfi'}; % magnet families for chromaticity matching
cLoptions.sext_fams     = {'sfi';'sfo';'sfm';'sdqd';'sdendq'}; % list of all sextupole families 
cLoptions.RBfams        = {'dqfa';'dqfb';'dqfc'}; % reverse bend families

cLoptions.eqfam = {'dip';'dip';'dip';'dip';'dip';'dip';'dip';'dip';'dip';...
                   'dip';'dip';'dip';'dip';'dip';'dip';'dip';'dip';'dip';...
                   'dipm';'dipm';'dipm';'dipm';'dipm';...
                   'Qf_Qfm';'Qf_Qfm';'Qf_Qfm';...
                   'dipm';'dipm';'dipm';'dipm';'dipm';'dipm';'dipm';...
                   'Sdend';'Sd';'Oxx_Oxy';'Oxx_Oxy';'Oyy';...
                   'Qfend_Qdend';'Qfend_Qdend';'Qf_Qfm';...
                   'Sfi_Sfo';'Sfm';'Sfi_Sfo'};
cLoptions.eqsca = ones(1,44);

load(strcat(erase(atroot,'atmat'),'/MAX4U/MagnetStrengthLimits.mat'));
load(strcat(erase(atroot,'atmat'),'/MAX4U/CandidateLattices/m4_standard/m4_standard.mat'));
ACHRO_ref = m4_standard.ACHROMAT;
V0 = 'auto';
bh =0.05;

%% Run cLatt options
m4UT = m4Uc_Latt(ACHRO,lattname,desc,cLoptions,ACHRO_ref,MagnetStrengthLimits,'V0',V0,'bh',bh,'basic');

% Below an example of how the 'plotLatt' function can be used to produce
% plots from the m4UT structure and save the results on a file.
%
% plotLatt(m4UT,'all','ymaxplot_dm',0.004,'zoom',2.0,'ymaxplot',0.004,'xminplot',-0.010,'xmaxplot',0.01,'dpminplotLMA',-0.25,'dpmaxplotLMA',0.25,'caxrange',[-10 0],'caxrange_r',[-10 -5],'nogrid');