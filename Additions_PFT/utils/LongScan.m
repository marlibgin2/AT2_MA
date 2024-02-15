function LongScan(i0,Nsteps,nfiles,MaxEmitt,MaxChro,MaxEta,choosef)
% Executes series of GALSS Scanns 
% The lattice, parameters to be varied and their ranges are defined in the
%   function glass_scan.
%
% inputs:
%    i0 : starting point - useful wen resuming from a previously interrupted scan
%    Nsteps: n. of steps in each variable
%    nfiles: number of sets of scans, each will be saved in one file
%    MaxEmit = maximum allowed emittance [pmrad]
%    MaxChrom = max allowed absolute value of chroaticities
%    Max Eta = max abslute value allowed dispersion at start of latice (centre of long
%               straight) [m]
%    choosef : if set to 'N', all stable lattices are recorded, regadless
%    of the limts set above. Otherwise thse limits (and Jx<3) are required
%    for a lattice to be reorded in the otu table
%
% Folder for saving data must exist before running this function. Folder must
% be in path. FOLDER chosen in a dialog
%
nvars = 7; % number of variables. CHANGE according to glass_scan.m
imax = Nsteps^nvars; % 
dI=round(imax/nfiles); % steps per file

folder=uigetdir;
  
k=1;
for i=i0:dI:imax
    datetime
    fprintf('Iteration n.: %12d  \n', k);
    fprintf('Starting lattice: %12d  \n', i);
    gscan1=glass_scan(Nsteps,dI,i,MaxEmitt,MaxChro,MaxEta,'Y',choosef,folder);
    k=k+1;
end
end

