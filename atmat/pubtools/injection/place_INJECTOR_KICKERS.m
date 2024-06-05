% ----------------------------------
% gets a RING or an ACHRomat, place 
% a MIK/DK in a chosen position
% expand the ACHRomat to a full RING
% equip with active cavities
% ----------------------------------
function RING = place_INJECTOR_KICKERS(Rin, DKposition,MIKposition)

% -------------------------------------
% is it an achromat? is it a full ring?
% -------------------------------------
RING = Rin;
spos = findspos(RING,1:length(RING)+1);
LRin = spos(end);
if LRin <=26.5
    RING = achromat2ring(RING);
end
spos = findspos(RING,1:length(RING)+1);

% -------------------------------
% possible positions for injector
% INJ1: 5.043m
% INJ2: 28.302m
% -------------------------------

%
% insert injector kicker 1
%

% -------------------------------------
% find nearby elements of INJ1 (5.043m)
% -------------------------------------
[INJ1_bmin INJ1_bmax] = INJborders(RING, DKposition);
% phc     = atmarker('PMH3','IdentityPass');
% pvc     = atmarker('PMV3','IdentityPass');
phc     = atcorrector('PMH3', 0.0, [0 0], 'CorrectorPass');
pvc     = atcorrector('PMV3',0.0, [0 0], 'CorrectorPass');

D       = spos(INJ1_bmax) - spos(INJ1_bmin);
FRAC1   = (DKposition-spos(INJ1_bmin))/D;
FRAC2   = (DKposition-spos(INJ1_bmin))/D;
RING = atinsertelems(RING,INJ1_bmin,FRAC1,phc,FRAC2,pvc);
spos = findspos(RING,1:length(RING)+1);
phc.Energy = 3.0e9;
phc.NumIntSteps=10;
phc.MaxOrder=3;
pvc.Energy = 3.0e9;
pvc.NumIntSteps=10;
pvc.MaxOrder=3;


% -------------------------------------
% find nearby elements of INJ1 (5.043m)
% -------------------------------------
[INJ2_bmin INJ2_bmax] = INJborders(RING, MIKposition);
[XMIK,   YMIK] = meshgrid(-0.012:0.001:0.012,0.006:-0.002:-0.006);
[PXMIK, PYMIK] = meshgrid(-0.012:0.001:0.012,0.006:-0.002:-0.006);
PXMIK = PXMIK.*0; PYMIK = PYMIK.*0; 

nn  = size(XMIK); nx=nn(2); ny=nn(1); 
mik = atepukick('MIK', nx, ny, XMIK, YMIK, PXMIK, PYMIK, 'ThinEPUPass');
mik.Energy = 3.0e9;
mik.NumIntSteps=10;
mik.MaxOrder=3;

D       = spos(INJ2_bmax) - spos(INJ2_bmin);
FRAC3   = (MIKposition-spos(INJ2_bmin))/D;
RING = atinsertelems(RING,INJ2_bmin,FRAC3,mik);

%
% introduce and activate cavities
%
RING = atenable_6d(RING);
RING = calculateGirderMaps(RING);
CAVi = findcells(RING,'FamName','CAV');
for ii=1:length(CAVi); RING{CAVi(ii)}.Voltage=1.1e6/length(CAVi); end  %%% 1.8e6/
for ii=1:length(CAVi); RING{CAVi(ii)}.TimeLag = 0.1843327465; end
hcmi = findcells(RING,'iscorH'); 
vcmi = findcells(RING,'iscorV');
sHCM = findspos(RING,hcmi);
sVCM = findspos(RING,vcmi);
end


function [bmin bmax] = INJborders(RING, sINJ)
% -------------------------------------
% find nearby elements of INJ1 (5.043m)
% -------------------------------------

%stest = 5.043;
spos = findspos(RING, 1:length(RING)+1);
[a bmax] =min(abs(spos-sINJ));
if spos(bmax)>sINJ
    smin=spos(bmax-1);
    bmin=bmax-1;
else
    smin = spos(bmax);
    bmin = bmax; 
    i = 0;
    while 1==1
        i=i+1;
        smax = spos(bmax+i);
        if smax>smin
            break
        end
    end
    bmax = bmax+i;
    bmin = bmax-1;
end
end
