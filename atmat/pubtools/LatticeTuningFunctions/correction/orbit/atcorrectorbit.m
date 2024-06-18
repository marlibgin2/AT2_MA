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
%     scalefactor,...   8) scale factor to correction (default: 0.75)
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



% % response matrix kicks
% kval=1e-5;
% delta=1e-3;

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

if nargin<7 || isempty(correctflags)
    correctflags=[true true];
end

if nargin<6 || isempty(neigSteerer)
%     neigSteerer=ones(10,1)*[length(indHCor)-20 length(indVCor)-20];
svhmax = min([numel(indHCor), numel(indBPM)]);
svvmax = min([numel(indVCor), numel(indBPM)]);
Niter = 10; fracStart=0.75;
scaleV = fracStart:(1-fracStart)/Niter:1;
svh = round(svhmax*scaleV);
svv = round(svvmax*scaleV);
neigSteerer = [svh; svv]';

% If corrector mean is being corrected using the RF, add another singular
% value to the horizontal plane
if correctflags(2)
    neigSteerer(:,1) = neigSteerer(:,1) + 1;
end

end

if nargin<8 || isempty(scalefactor)
    if printouttext, disp(' --- scale set to .75'); end
    scalefactor=0.75;
end

if nargin<9 || isempty(ModelRM)
    if printouttext, disp(' --- computing Orbit Response Matrix (ORM)'); end
    ModelRM=[];
end

if nargin<10 || isempty(reforbit)
    if printouttext, disp(' --- reference orbit = 0'); end
    reforbit=zeros(2,length(indBPM));
end

if scalefactor<0 || scalefactor>1
    if printouttext
        disp(' --- scale factor out of range. Set to .75'); end
    scalefactor=0.75;
end

% if correctflags(1) % dpp correction
%     rmsel=[1 2 3];
% else
%     rmsel=[1 2];
% end



% Check whether lattice is 6D or 4D
use6d = check_6d(rerr);

% Determine corrector model type and the correct field to use
% NB! The assumption is that this does not vary across corrector families
if strcmpi(rerr{indHCor(1)}.PassMethod,'CorrectorPass'), xfname = 'KickAngle'; xfi = 1; xfj = 1; ...
else, xfname = 'PolynomB'; xfi = 1; xfj = 1; ...
end
if strcmpi(rerr{indVCor(1)}.PassMethod,'CorrectorPass'), yfname = 'KickAngle'; yfi = 1; yfj = 2; ...
else, yfname = 'PolynomA'; yfi = 1; yfj = 1; ...
end



% Attempt to get initial orbit
% NB! For some error seeds no initial orbit exists. As we will attempt to
% deal with that below, the warning of an ill-conditioned matrix is
% temporarily disabled.
warning('OFF','MATLAB:illConditionedMatrix');
if use6d
    o=findorbit6Err(rerr,indBPM,inCOD);
else
    o=findorbit4Err(rerr,0,indBPM,inCOD);
end
warning('ON','MATLAB:illConditionedMatrix');


% NB! Before continuing a check is needed to avoid getting into a recursive
% "lock". Throw an error if the recursive depth is exceeded.
maxRecursiveDepth = 2;
if getRecursiveDepth > maxRecursiveDepth
    error('atcorrectorbit:No closed orbit could be found!');
end

% If no orbit was found, a good guess of initial corrector settings is
% needed.
errorFraction = 0.1;
if any(isnan(o(:)))
    if printouttext, fprintf('No orbit found! Searching for initial corrector setting ...\n'); end
    RING_RED = rerr;

    % Zero the correctors for the reduced lattice
    RING_RED = atsetfieldvalues(RING_RED,indHCor,xfname,{xfi, xfj},0);
    RING_RED = atsetfieldvalues(RING_RED,indVCor,yfname,{yfi, yfj},0);

    % Reduce the translation errors by a significant factor (errorFraction)
    % As a large fraction of the orbit errors are usually due to quad
    % misalignments, hopefully this allows the calculation of a closed
    % orbit. 
    I = atgetcells(rerr,'T1');
    redVals.T1 = cellfun(@(x) x*errorFraction,atgetfieldvalues(rerr(I),'T1'),'UniformOutput',false);
    redVals.T2 = cellfun(@(x) x*errorFraction,atgetfieldvalues(rerr(I),'T2'),'UniformOutput',false);
    RING_RED = atsetfieldvalues(RING_RED,I,'T1',redVals.T1);
    RING_RED = atsetfieldvalues(RING_RED,I,'T2',redVals.T2);

    % Attempt to correct the easier case
    % NB! This is a recursive call!
    RING_RED = atcorrectorbit(RING_RED,indBPM,indHCor,indVCor,[],[],[],[],ModelRM);

    % Extract the obtained corrections and scale them up, then use them as
    % the initial guess for the original lattice.
    vals.hcm = atgetfieldvalues(RING_RED,indHCor,xfname,{xfi, xfj});
    vals.vcm = atgetfieldvalues(RING_RED,indVCor,yfname,{yfi, yfj});
    rerr = atsetfieldvalues(rerr,indHCor,xfname,{xfi, xfj},vals.hcm/errorFraction);
    rerr = atsetfieldvalues(rerr,indVCor,yfname,{yfi, yfj},vals.vcm/errorFraction);

    % Re-calculate the closed orbit at the BPMs
    if use6d
        o=findorbit6Err(rerr,indBPM,inCOD);
    else
        o=findorbit4Err(rerr,0,indBPM,inCOD);
    end
end

ox0=o(1,:);
oy0=o(3,:);

% Get BPM weight information.
W = cell2mat(atgetfieldvalues(rerr,indBPM,'Weight'));
W(isnan(W)) = 1;    % Any BPMs without specified weight are assumed to have weight 1.

% Compute momentum compaction
alpha=mcf(rerr);

% Load or compute response matrix
% NB! This is done outside the orbit correction loop for speed reasons; the
% loss of accuracy is generally not an issue as the correction loop can
% tolerate fairly large errors in the RM. 
[RMH, RMV] = calcRM(ModelRM);

% if isempty(ModelRM)
%     % get orbit RM
%     if printouttext
%         disp('Calculating ORM using adapted LOCO-routine'); end
%         ModelRM = getlinearrespmat(rerr,indBPM,indHCor,indVCor);
% 
%     %     ModelRM=getresponsematrices(...
%     %         rerr,...          %1 AT lattice
%     %         indBPM,...      %2 bpm indexes in at lattice
%     %         indHCor,...     %3 h cor indexes
%     %         indVCor,...     %4 v cor indexes
%     %         [],...     %5 skew cor indexes
%     %         [],...     %6 quad cor indexes
%     %         [],...     %7 sext cor indexes
%     %         inCOD,...       %8 initial coordinates
%     %         rmsel...      %9 specifiy rm to be computed
%     %         );
% end


    function [RMH, RMV] = calcRM(ModelRM)

        if isempty(ModelRM)
            % get orbit RM
            if printouttext
                disp('  Calculating ORM using adapted LOCO-routine'); end
            ModelRM = getlinearrespmat(rerr,indBPM,indHCor,indVCor);
        end

        if ~correctflags(1) % dpp correction
            ModelRM.OrbHDPP=[];
            ModelRM.OrbVDPP=[];
        end

        ormH=ModelRM.OrbHCor;
        ormV=ModelRM.OrbVCor;


        % Re-format the calculated ORM for the orbit correction calculations
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

        % Rescale the ORM to take BPM weights into account
        for n = 1:numel(indBPM)
            RMH(n,:) = RMH(n,:).*W(n,1);
            RMV(n,:) = RMV(n,:).*W(n,2);
        end

        % Get the SV sizes
        [Uh, Sh, Vh] = svd(RMH);
        [Uv, Sv, Vv] = svd(RMV);
        sh = diag(Sh);
        sv = diag(Sv);

%         % Then truncate the very small singular values, as these represent
%         % a mode that cannot be corrected
%         SVthreshold = 1e-8;
%         shm = find(sh./max(sh) > SVthreshold,1,'last');
%         svm = find(sv./max(sv) > SVthreshold,1,'last');
%         neigSteerer(neigSteerer(:,1) > shm,1) = shm;
%         neigSteerer(neigSteerer(:,2) > svm,1) = svm;

        % Then apply RM regularization, i.e. rescale the SVs for better
        % numerical behaviour
        hlambda = max(sh) * 1e-3;
        vlambda = max(sv) * 1e-3;
        shn = (sh + hlambda).^2 ./ sh; for n = 1:numel(shn), Sh(n,n) = shn(n); end
        svn = (sv + vlambda).^2 ./ sv; for n = 1:numel(svn), Sv(n,n) = svn(n); end
        RMH = Uh*Sh*Vh';
        RMV = Uv*Sv*Vv';
    end


%% MAIN CORRECTION LOOP
orbitThreshold = 1e-6;
convergenceThreshold = orbitThreshold / 2;
iter = 0;

% Get starting orbit at entrance (for later tune computation and guess
% for next orbit search) and all BPMs
if use6d
    o=findorbit6Err(rerr,[1; indBPM(:)],inCOD);
else
    o=findorbit4Err(rerr,0,[1; indBPM(:)],inCOD);
end

% Loop logic:
% a) Ramp up 
% b) Then, if rms orbit above threshold, continue with last set of SVs until:
%     1. rms orbit below threshold (i.e. orbit is good enough)
%     2. rms orbit improvement below threshold (i.e. not worth continuing)
while true  %iter=1:Niter
    iter = iter + 1;
    if printouttext
        disp(['Orbit correction iter ' num2str(iter,'%d, ') 'n-eig: ' num2str(neigSteerer(min([iter, size(neigSteerer,1)]),:),'%d, ')]);
    end

    % Store previous orbit and update best corrector settings
    if iter > 1
        oxprev = ox;
        oyprev = oy;
    end

    % Get current corrector strengths
    corh0=atgetfieldvalues(rerr,indHCor,xfname,{xfi,xfj});
    corv0=atgetfieldvalues(rerr,indVCor,yfname,{yfi,yfj});

    % Get current orbit
    ox=o(1,2:end);
    oy=o(3,2:end);

    % Subtract reference orbit
    ox=ox-reforbit(1,:);
    oy=oy-reforbit(2,:);

    if printouttext
        disp(['  Orbit RMS value:  ' num2str(std(ox)*1e6,'%.2f') ' / ' num2str(std(oy)*1e6,'%.2f µm')]);
    end
    if printouttext && iter > 1
        disp(['  Orbit RMS change: ' num2str((std(ox)-std(oxprev))*1e6,'%.2f µm') ' / ' num2str((std(oy)-std(oyprev))*1e6,'%.2f µm')]);
    end

    % If threshold requirement is met, exit the loop
    if std(ox) <= orbitThreshold && std(oy) <= orbitThreshold
        break;
    end

    % Temporary increases in orbit error when increasing the SVs
    % are to be expected, so don't exit during the SV ramp plus a few
    % iterations
    if iter > size(neigSteerer,1) + 2
        % If not enough improvement in rms orbit, exit loop
        if std(oxprev) - std(ox) < convergenceThreshold && ...
                std(oyprev) - std(oy) < convergenceThreshold
            break;
        end
    end

    % Apply BPM weights
    oxw = ox.*W(:,1)';
    oyw = oy.*W(:,2)';


    % Compute correction
    if correctflags(2) % mean 0
        dch=qemsvd_mod(RMH,[oxw';0],neigSteerer(min([iter, size(neigSteerer,1)]),1));
        dcv=qemsvd_mod(RMV,[oyw';0],neigSteerer(min([iter, size(neigSteerer,1)]),2));
    else % no constraint on correctors mean
        dch=qemsvd_mod(RMH,oxw',neigSteerer(min([iter, size(neigSteerer,1)]),1));
        dcv=qemsvd_mod(RMV,oyw',neigSteerer(min([iter, size(neigSteerer,1)]),2));
    end


    % Get total correctors values and apply scaling
    if correctflags(1)
        hs=corh0-dch(1:end-1)*scalefactor;
        vs=corv0-dcv(1:end-1)*scalefactor;
        % energy deviation
        dd=-dch(end)*scalefactor;
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
            disp(['  Delta RF: ' num2str(-alpha*(dd)*f0) ' Hz']);
        end
    end

    % Make an educated guess for the new closed orbit
    % Comment: will likely not work as well with a non-zero reference
    inCOD(1:4) = inCOD(1:4)*(1-scalefactor);
    
    % Get the orbit at entrance (for later tune computation and guess
    % for next orbit search) and all BPMs 
    if use6d
        o=findorbit6Err(rtest,[1; indBPM(:)],inCOD);
    else
        o=findorbit4Err(rtest,0,[1; indBPM(:)],inCOD);
    end

    % Update the input now that we know it
    inCOD = o(:,1);

    % Attempt to calculate the tunes to check lattice stability
    try
        t=tunechrom(rtest,'orbit',inCOD);   % Re-use the calculated orbit
    catch
        t = nan(2,1);
    end
    if printouttext
        disp(['  [nux, nuy]: ' num2str(t(1)) ' / ' num2str(t(2)) ]);
    end
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

    % Update the response matrix, if not given
    % NB! Expensive calculation-wise but may be required if the orbit has
    % changed a lot.
    if isempty(ModelRM)
        [RMH, RMV] = calcRM(ModelRM);
    end

end

%% DATA EXTRACTION AND POST-PROCESSING

% Get current orbit
if use6d
    o=findorbit6Err(rcor,indBPM,inCOD);
    dpp = o(5,1);
else
    o=findorbit4Err(rcor,0,indBPM,inCOD);
    dpp = 0;
end
oxc=o(1,:);
oyc=o(3,:);

% Get current corrector strengths
Lh=atgetfieldvalues(rcor,indHCor,'Length');
Lv=atgetfieldvalues(rcor,indVCor,'Length');
Lh(Lh == 0) = 1;    % In case of thin correctors, KickAngle should already have been used
Lv(Lv == 0) = 1;    % In case of thin correctors, KickAngle should already have been used

hsL=hs.*Lh;
vsL=vs.*Lv;

if printouttext
    % display results
    disp(' ');
    disp(' ');
    disp(['      before' '    ' '-->' '    ' 'after'])
    disp(['oX: ' num2str(std(ox0-reforbit(1,:))*1e6,'%3.3f') ' -> ' num2str(std(oxc-reforbit(1,:))*1e6,'%3.3f') 'um']);
    disp(['oY: ' num2str(std(oy0-reforbit(2,:))*1e6,'%3.3f') ' -> ' num2str(std(oyc-reforbit(2,:))*1e6,'%3.3f') 'um']);
    disp(['    ' 'min' '    ' 'mean' '    ' 'max'])
    disp(['hs:'  num2str([min(hsL) mean(hsL) max(hsL)]*1e3,' %2.2f ') ' mrad'])
    disp(['vs:'  num2str([min(vsL) mean(vsL) max(vsL)]*1e3,' %2.2f ') ' mrad'])
    disp(['dpp @ s=0: ' num2str(dpp)])
end
end


%% -----HELPER FUNCTIONS-----
function N = getRecursiveDepth
stack = dbstack;
N = sum(strcmp({stack.name},stack(2).name));
end