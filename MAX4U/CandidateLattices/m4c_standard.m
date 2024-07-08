%% m4c_standard
% This script is a template to set up workspace variables needed
% to run the "m4_cLatt" script, whuch in turn, runs all "cLatt" function 
% options in a series of steps creating/overwriting a structure named
% "m4T" in the matlab workspace, saving this variable
% in a "m4T.mat" file and saving all matlab output in a logfile
%
% The intended workflow is 
%  1. Edit the statements in the "Lattice Specific Data" section below.
%  2. Run the script.
%  3. Rename or copy the resulting "m4T" sctructure accordingly.
%  4. Execution may be interrupted at any moment and the most recent 
%     m4T scruture can be recovered by loading it from the file 
%     saved on disk. From that point on, further calculation and 
%     updates to the m4T structure may be carried out by running 
%     cLatt with specific input options. See the m4_cLatt script for 
%     many possible calls to cLatt.
%  
% The structure "cLoptions" is stored as a field in the "m4T" structure 
% and contains various lattice evaluation optional settings. These are
% listed in the "General initialisation" section of the m4_cLatt script 
% (which sets their default values). If values different from the defauls 
% are desired for fields in the cLoptions structure, these can be changed 
% at anytime, e.g., by directly editing the correspondig fields or by 
% typing:
%    >>cLpptions = m4T.cLoptions; % extracts cLoptions to the workspace.
%    >>cLoptions.DAoptions.nturns = 2028 % implements the desired changes -
%    in this example, we change the number of turns used in dynamic
%    aperture calcualations
%
%    >>m4T=cLatt(m4T, 'cLoptions',cLoptions,'verbose',1); % stores the
%                                          changes in the "m4T" structure.
%    Note that such changes do not automatically re-run any evaluations 
%    that are affected by them - these need to be rerun explicitly 
%    to update the corresponding fields in the m4T structure, e.g.
%
%    >> m4T = cLatt(m4T,'DAxy','verbose',2);
%    
%
%% Lattice specific data
lattname = 'm4 Standard Lattice';
desc = 'm4-20230102_fromtracy20121107-430-bare.opa';
diary_file=strcat(lattname,'_log_',datestr(now,30));
diary(diary_file);

m4U('a1','Full','rbk',0.0);
ACHRO = ACHROGRD_a1; % The cell array with the AT2 lattice to be evaluated

clear('cLoptions');
cLoptions.All_famsO=LatticeOptData.All_famsO; %   optional, If empty m4_cLatt finds out the magnet 
%                           family names. Use this in case a specific order
%                           of family names is desired.
cLoptions.ringtune_fams = LatticeOptData.ringtune_fams;   % magnet families for tune matching
cLoptions.chrom_fams    = LatticeOptData.chrom_fams; % magnet families for chromaticity matching
cLoptions.eqfam         = LatticeOptData.eqfam;

cLoptions.eqsca         = LatticeOptData.eqsca;
cLoptions.ErrorModel    = errormodel_DDRchallenging('gdran',1.0,...
                            'mgalran',1.0,'mulsys',1.0,'mulran',1.0, ...
                            'strran',1.0,'bpmran',1.0);


%% Run cLatt options
m4_cLatt;

diary off