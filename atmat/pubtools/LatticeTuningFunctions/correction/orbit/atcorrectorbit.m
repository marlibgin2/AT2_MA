function [rcor,inCOD,hs,vs]=atcorrectorbit(...
    rerr,...
    indBPM,...
    indHCor,...
    indVCor,...
    inCOD,...
    neigSteerer,...
    correctflags,...
    scalefactor,...
    ModelRM,...
    reforbit,...
    steererlimit,...
    printouttext)
%
% Closed orbit correction.
%
% function [...
%    rcor,...           1) corrected lattice
%    inCOD,...          2) initial COD (dpp is stored here)
%    hs,vs...           3,4) total steerers strengths after correction
%    ]=atcorrectorbit(...
%     rerr,...          1) AT lattice to correct
%     indBPM,...        2) Nbx1 bpm indexes
%     indHCor,...       3) Nhx1 hor. cor indexes
%     indVCor,...       4) Nvx1 ver. cor indexes
%     inCOD,...         5) 6x1 initial COD guess
%     neigSteerer,...   6) 2xNiter eigenvectors for correction H and V at
%                          each iteration (default: [Nh/2 Nv/2])
%     correctflags,...  7) correct [dpp mean0](default: [true true])
%     scalefactor,...   8) scale factor to correction (default: 1.0)
%     ModelRM,...       9) ModelRM.Orb(H/V)Cor = 4x1 cell of orbit response matrix
%                          ModelRM.Orb(H/V)DPP = 6x1 array of orbit
%                          response to dpp
%                          if [] compute RM (default: [])
%     reforbit,...      10) 2xNbpm reference orbit to correct to (default 0*2xNb)
%     steererlimit      11) 2x1 limit of steerers abs(steerer)<steererlimit
%                           (default: [], no limits)
%     printouttext      12) if 1 or true, display rms orbit
%     )
%
% features impelemented:
% limit correctors strengths
% ddp correction
% sum of steerers = 0
% 6D orbit with BPM errors
% initial coordinate
% correction to reference orbit refx refy
% use atsetfieldvalues, atgetcells
%
%
%see also: qemsvd_mod findorbit6Err getresponsematrices



% response matrix kicks
kval=1e-5;
delta=1e-3;

alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
f0=rerr{indrfc(1)}.Frequency;

% default arguments
if nargin<12 || isempty(printouttext)
    printouttext=true;
end
if nargin<11 || isempty(steererlimit)
    steererlimit=[];
end

if nargin<2 || isempty(indBPM)
    if printouttext
        disp('No BPM indices, guessing based on lattice'); end
    indBPM=find(atgetcells(rerr,'Class','Monitor'));
end

if nargin<3 || isempty(indHCor)
    if printouttext
        disp('No hor. corrector indices, guessing based on lattice'); end
    indHCor=find(atgetcells(rerr,'iscorH','H'));
end

if nargin<4 || isempty(indVCor)
    if printouttext
        disp('No ver. corrector indices, guessing based on lattice'); end
    indVCor=find(atgetcells(rerr,'iscorV','V'));
end


if nargin<5 || isempty(inCOD)
    inCOD=[0 0 0 0 0 0]';
end

if nargin<6 || isempty(neigSteerer)
    neigSteerer=ones(10,1)*[length(indHCor)-20 length(indVCor)-20];
end

if nargin<7 || isempty(correctflags)
    correctflags=[true true];
end

if nargin<8 || isempty(scalefactor)
    if printouttext
        disp(' --- scale set to 1.0'); end
    scalefactor=0.2;
end

if nargin<9 || isempty(ModelRM)
    if printouttext, disp(' --- computing orbit Response matrix'); end
    ModelRM=[];
end

if nargin<10 || isempty(reforbit)
    if printouttext, disp(' --- reference orbit = 0'); end
    reforbit=zeros(2,length(indBPM));
end

if scalefactor<0 || scalefactor>1
    if printouttext
        disp(' --- scale factor out of range. Set to 1.0'); end
    scalefactor=0.2;
end

if correctflags(1) % dpp correction
    rmsel=[1 2 3];
else
    rmsel=[1 2];
end

% load or compute response matrix
if isempty(ModelRM)
    % get orbit RM
    if printouttext
        disp('get orbit RM'); end
    
%     ModelRM=getresponsematrices(...
%         rerr,...          %1 AT lattice
%         indBPM,...      %2 bpm indexes in at lattice
%         indHCor,...     %3 h cor indexes
%         indVCor,...     %4 v cor indexes
%         [],...     %5 skew cor indexes
%         [],...     %6 quad cor indexes
%         [],...     %7 sext cor indexes
%         inCOD,...       %8 initial coordinates
%         rmsel...      %9 specifiy rm to be computed
%         );
ModelRM = getlinearrespmat(rerr,indBPM,indHCor,indVCor);    
    
    if ~correctflags(1) % dpp correction
        ModelRM.OrbHDPP=[];
        ModelRM.OrbVDPP=[];
    end
    
end

ormH=ModelRM.OrbHCor;
ormV=ModelRM.OrbVCor;
% kval=ModelRM.kval;
% dppH=ModelRM.OrbHDPP;
% dppV=ModelRM.OrbVDPP;
% delta=ModelRM.delta;

% get initial orbit
o=findorbit6Err(rerr,indBPM,inCOD);
ox0=o(1,:);
oy0=o(3,:);

%rerr0=rerr;

% Determine corrector model type and the correct field to use
if strcmpi(rerr{indHCor(1)}.PassMethod,'CorrectorPass'), xfname = 'KickAngle'; xfi = 1; xfj = 1; ...
else, xfname = 'PolynomB'; xfi = 1; xfj = 1; ...
end
if strcmpi(rerr{indVCor(1)}.PassMethod,'CorrectorPass'), yfname = 'KickAngle'; yfi = 1; yfj = 2; ...
else, yfname = 'PolynomB'; yfi = 1; yfj = 1; ...
end

% Get BPM weight information.
W = cell2mat(atgetfieldvalues(rerr,indBPM,'Weight'));
W(isnan(W)) = 1;    % Any BPMs without specified weight are assumed to have weight 1.

% iterate correction
Niter=size(neigSteerer,1);
for iter=1:Niter
    
    if printouttext
        disp(['Orbit correction iter ' num2str(iter,'%d, ') 'n-eig: ' num2str(neigSteerer(iter,:),'%d, ')]);
    end
    
    % initial corrector strengths
    corh0=atgetfieldvalues(rerr,indHCor,xfname,{xfi,xfj});
    corv0=atgetfieldvalues(rerr,indVCor,yfname,{yfi,yfj});
    
    % get current orbit
    o=findorbit6Err(rerr,indBPM,inCOD);
    ox=o(1,:);
    oy=o(3,:);
    
    % subtract reference orbit
    ox=ox-reforbit(1,:);
    oy=oy-reforbit(2,:);
    



    % build RMs
    if correctflags(1) && correctflags(2) % dpp and mean0
        dppH=ModelRM.OrbHDPP;
        dppV=ModelRM.OrbVDPP;
        RMH=[ [ormH{1};ones(size(indHCor(:)))'] [dppH(:);0] ];
        RMV=[ [ormV{3};ones(size(indVCor(:)))'] [dppV(:);0] ];
    elseif correctflags(1) && ~correctflags(2)% dpp no mean 0
        dppH=ModelRM.OrbHDPP;
        dppV=ModelRM.OrbVDPP;
        RMH=[ ormH{1} dppH(:) ];
        RMV=[ ormV{3} dppV(:) ];
    elseif ~correctflags(1) && correctflags(2) % mean0 no dpp
        RMH=[ormH{1};ones(size(indHCor(:)))'];
        RMV=[ormV{3};ones(size(indVCor(:)))'];
    elseif ~correctflags(1) && ~correctflags(2) % no dpp no mean0
        RMH=ormH{1};
        RMV=ormV{3};
    end
    
    % Apply BPM weights
    ox = ox.*W(:,1)';
    oy = oy.*W(:,2)';

    % Rescale the ORM
    for n = 1:numel(indBPM)
        RMH(n,:) = RMH(n,:).*W(n,1);
        RMV(n,:) = RMV(n,:).*W(n,2);
    end
    
    % compute correction
    if correctflags(2) % mean 0
        dch=qemsvd_mod(RMH,[ox';0],neigSteerer(iter,1));
        dcv=qemsvd_mod(RMV,[oy';0],neigSteerer(iter,2));
    else % no constraint on correctors mean
        dch=qemsvd_mod(RMH,ox',neigSteerer(iter,1));
        dcv=qemsvd_mod(RMV,oy',neigSteerer(iter,2));
    end
    
    
    % get total correctors values and apply scaling
    if correctflags(1)
        hs=corh0-dch(1:end-1)*scalefactor;
        vs=corv0-dcv(1:end-1)*scalefactor;
        % energy deviation
        dd=-dch(end);
    else
        hs=corh0-dch*scalefactor;
        vs=corv0-dcv*scalefactor;
    end
    
    % limit steerers strengths
    if ~isempty(steererlimit)
        hs(abs(hs)>steererlimit(1))=steererlimit(1);
        vs(abs(vs)>steererlimit(2))=steererlimit(2);
    end
    
    rtest=atsetfieldvalues(rerr,indHCor,xfname,{xfi,xfj},hs);
    rtest=atsetfieldvalues(rtest,indVCor,yfname,{yfi,yfj},vs);
    
   
    if correctflags(1)
        
        rtest=atsetfieldvalues(rtest,indrfc,'Frequency',f0-alpha*(dd)*f0);
        
        if printouttext
            disp(['Delta RF : ' num2str(-alpha*(dd)*f0) ' Hz']);
        end
    end
    
    %[~,t,~]=atlinopt(rtest,0,1);
    t=tunechrom(rtest,0);
    if not(isnan(t(1)) || isnan(t(2)))
        % apply correction in lattice
        rcor = rtest;
    else 
        rcor = rerr;
        if printouttext
            disp('Orbit correction fails, stop at previous step.')
        end
        break
    end
    
    % lattice start point for next iteration
    rerr=rcor;
end

% get current orbit
o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);
    oxc = oxc.*W(:,1)';
    oyc = oyc.*W(:,2)';

Lh=atgetfieldvalues(rcor,indHCor,'Length');
Lv=atgetfieldvalues(rcor,indVCor,'Length');
hsL=hs.*Lh;
vsL=vs.*Lv;

if printouttext
    % display results
    disp(['      before' '    ' '-->' '    ' 'after'])
    disp(['oX: ' num2str(std(ox0-reforbit(1,:))*1e6,'%3.3f') ' -> ' num2str(std(oxc-reforbit(1,:))*1e6,'%3.3f') 'um']);
    disp(['oY: ' num2str(std(oy0-reforbit(2,:))*1e6,'%3.3f') ' -> ' num2str(std(oyc-reforbit(2,:))*1e6,'%3.3f') 'um']);
    disp(['    ' 'min' '    ' 'mean' '    ' 'max'])
    disp(['hs:'  num2str([min(hsL) mean(hsL) max(hsL)]*1e3,' %2.2f ') ' mrad'])
    disp(['vs:'  num2str([min(vsL) mean(vsL) max(vsL)]*1e3,' %2.2f ') ' mrad'])
    disp(['dpp: ' num2str(inCOD(5))])
end
