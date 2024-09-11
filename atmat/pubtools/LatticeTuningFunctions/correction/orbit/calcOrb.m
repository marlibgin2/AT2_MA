function [RINGc,orb0,orb,hcor,vcor] = calcOrb(varargin)
% Calculates, plots and corrects the closed orbit
% 
% This is a higher level wrapper function
% 
%% Inputs
% Mandatory input arguments
% RING : AT2 lattice array
%
% Optional input arguments
% ORM                   : orbit reponse matrix
% verbose               : defines level of verbose output, default=0, i.e. no output
% OCoptions             : structure wiht the fields
% OCoptions.inCOD       : inital guess for the orbit
% OCoptions.neigen      : 2xNiter eigenvectors for correction H and V at
%                         each iteration (default: [Nh/2 Nv/2])
% OCoptions.cflags      : correct [dpp mean0](default: [true true])
% OCoptions.scale       : scale factor to correction (default: 0.75)
% OCoptions.reforbit    : 2xNbpm reference orbit to correct to (default 0*2xNb)
% OCoptions.steererlimit       : 2x1 limit of steerers abs(steerer)<steererlimit
%                           (default: [], no limits)
% Optional flags
% correct: corrects the orbit
% plot : plots the orbit
%
%% Outputs
% RINGc: corrected ring (if correction is not asked for, this is the same as the input lattice)
% orb0: (nx2) array: (X,Y) Initial orbit  [m] 
% orb: (nx2) array: (X,Y) Corrected orbit [m] (only if correction is done)
% 
%% Usage examples
% [RINGc, orb0, orb] = calcOrb(RING,'plot','correct');
% [~, orb0, orb]     = calcOrb(RING,'plot');
% RINGc = calcOrb(RING,'plot','correct');
% calcOrb(RING,'plot');

%% History
% PFT 2024/03/02, first version
% PFT 2024/03/10: Bug fix to get the correct output lattice
% PFT 2024/06/30: added search for corrector indices based on family names
% PFT 2024/07/08: added search for BPM indices based on alternative family
%                 names
% PFT 2024/07/16: improved handling of verbose option
% PFT 2024/07/25: added possibility of fixing the orbit reponse matrix
% PFT 2024/07/31: improved handling of correctf=false
%
%% Input argument parsing
RING           = getargs(varargin,[]);
ORM            = getoption(varargin,'ORM',[]);
plotf          = any(strcmpi(varargin,'plot'));
correctf       = any(strcmpi(varargin,'correct'));
verboselevel   = getoption(varargin,'verbose',0);
OCoptions      = getoption(varargin,'OCoptions',struct());

if (isempty(fields(OCoptions)))
    OCoptions.inCOD          = [];
    OCoptions.neigen         = [];
    OCoptions.cflags         = [];
    OCoptions.scale          = 0.75;
    OCoptions.reforbit       = [];
    OCoptions.steererlimit   = [];
end

inCOD        = OCoptions.inCOD;
neigen       = OCoptions.neigen;
cflags       = OCoptions.cflags;
scale        = OCoptions.scale;
reforbit     = OCoptions.reforbit;
steererlimit = OCoptions.steererlimit;

%% Calculates the closed orbit
% View the orbit, including BPM errors
setoption('WarningDp6D',false); % avoids warning messages
iBPM = findcells(RING,'FamName','BPM');
if (isempty(iBPM))
    iBPM=findcells(RING,'FamName','mon');
end
RINGc = RING;
if (isempty(iBPM))
    fprintf('%s Error in calcOrb; no BPMs in Lattice to plot orbit, aborting... \n', datetime);
    orb0=nan;
    orb=nan;
    return
end
sBPM = findspos(RING,iBPM);
orb0 = findorbit6Err(RING,iBPM);

if (plotf)
    figure; plot(sBPM,1e3*orb0([1 3],:)); xlim([0 528]);
            xlabel('s [m]'); ylabel('x,y [mm]');grid;legend('X','Y');
            title('Before correction');
end

%% Corrects the orbit
if (correctf)
    indHCor=find(atgetcells(RING,'iscorH','H'));
    if (isempty(indHCor))
        indHCor=findcells(RING,'FamName','ch');
    end
    indVCor=find(atgetcells(RING,'iscorV','V'));
    if (isempty(indVCor))
        indVCor=findcells(RING,'FamName','cv');
    end
    sHcor= findspos(RING,indHCor);
    sVcor= findspos(RING,indVCor);

%    RINGc = atcorrectorbit(RING,[],[],[],[],[140 120; 160 140; 180 160; ...
%                           ones(10,1)*[200 180]],[true true],0.75,...
%                           [],[],[0.38, 0.38]*1e-3,verbosef);
     
     [RINGc,~,hcor,vcor]= atcorrectorbit(RING,iBPM,indHCor,indVCor,inCOD,neigen,...
                           cflags,scale,ORM,reforbit,steererlimit,(verboselevel-1)>0);
     
    % Calculate the new orbit and plot in the former figure
    orb = findorbit6Err(RINGc,iBPM);
    if (plotf)
        figure;plot(sBPM,1e6*orb([1 3],:));xlim([0 528]); 
               xlabel('s [m]'); ylabel('x,y [µm]');grid;legend('X','Y');
               title('After correction');

        figure;plot(sHcor,hcor*1000);hold on; plot(sVcor,vcor*1000);
        legend('Hcor','Vcor');xlabel('s [m]'); ylabel('Cor[mrad]');
        grid;title('Corrector Strengths')

    end   
else
    orb=orb0;
end





