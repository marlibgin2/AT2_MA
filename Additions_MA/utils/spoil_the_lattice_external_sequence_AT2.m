function [RINGs, RING, MAGe, GIRe, Xbpm, Ybpm] = spoil_the_lattice_external_sequence_AT2(RING,nA)%, OrbitCorrectionFlag)
T = load('ModelRM.mat'); ModelRM = T.ModelRM; clear T
if nargin<2
    nA = 1; % number of achromats
end
% -------------------------
% spoil the entire lattice, 
% -------------------------
% assign zero shifts to all elements, creating the T1, T2 fields
allelemi = 1: length(RING);
RING     = atsetshift(RING,allelemi,0,0);
RING     = atsettilt(RING,allelemi,0);
spos     = findspos(RING,allelemi);

% ----------------------------------------
% single magnet error table (RMS) --> MAGe
% ----------------------------------------
%      grad(frac)   dx(um)    dy(um)
gradZero  = 1;  % 1 turn off/on the gradient errors
shiftZero = 1e-1;  % 1 turn off/on the displacement errors
rollZero  = 0;  % 1 turn off/on the roll errors
eQ = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % quadrupole
eR = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % reverse-bends
eS = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % sextupole
eO = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % octupole
eD = [ 5e-4*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % dipole
MAGe.eQ = eQ; MAGe.eR = eR; MAGe.eS = eS; MAGe.eO = eO; MAGe.eD = eD;

% -------------------------------
% girder random error table (RMS)
% -------------------------------
%        sway(um) heave(um) yaw(urad) pitch(urad) roll(urad)
grdZero = 0;
eGr{1}  = [100      100       25        25          10] * 1e-6  *grdZero; 
eGr{2}  = [100      100       25        25          10] * 1e-6  *grdZero; 
eGr{3}  = [100      100       25        25          10] * 1e-6  *grdZero; 
eGr{4}  = [100      100       25        25          10] * 1e-6  *grdZero; 
eGr{5}  = [100      100       25        25          10] * 1e-6  *grdZero; 
eGr{6}  = [100      100       25        25          10] * 1e-6  *grdZero; 
eGr{7}  = [100      100       25        25          10] * 1e-6  *grdZero; 

% --------------------------------------
% girder test error table - fixed values
% --------------------------------------
%        sway(um) heave(um) yaw(urad) pitch(urad) roll(urad)
egt{1}  = [  0        0        0         0          0] * 1e-6; 
egt{2}  = [300        0      100         0          0] * 1e-6; 
egt{3}  = [123        0     -100         0          0] * 1e-6; 
egt{4}  = [  0        0        0         0          0] * 1e-6; 
egt{5}  = [  0     -123        0        50          0] * 1e-6; 
egt{6}  = [  0     -300        0       -50          0] * 1e-6; 
egt{7}  = [  0        0        0         0          0] * 1e-6; 

girder_move_type = 'nomove'; % 'fixed' / 'nomove' / 'random'
GIRe.type = girder_move_type;
if strcmpi(girder_move_type,'random')   
    GIRe.gir = eGr;
elseif strcmpi(girder_move_type,'fixed') 
    GIRe.gir = egt;
else
    GIRe.gir = [];
end

% ---------------
% AT quad indices
% ---------------
q1i  = findcells(RING,'FamName','Q1');
q2i  = findcells(RING,'FamName','Q2');
q3i  = findcells(RING,'FamName','Q3');
q4i  = findcells(RING,'FamName','Q4');
q5i  = findcells(RING,'FamName','Q5');
qi   = sort([q1i q2i q3i q4i q5i]); qcol=[255 100 0]/255;
r1i  = findcells(RING,'FamName','R1');
r2i  = findcells(RING,'FamName','R2');
ri   = sort([r1i r2i]);rcol=[255 0 200]/255;

% ---------------
% AT sext indices
% ---------------
s1i    = findcells(RING,'FamName','S1');
s2i    = findcells(RING,'FamName','S2');
s3i    = findcells(RING,'FamName','S3');
s4i    = findcells(RING,'FamName','S4');
s5i    = findcells(RING,'FamName','S5');
si     = sort([s1i s2i s3i s4i s5i]); scol=[0 255 200]/255;

% ---------------
% AT oct indices
% ---------------
o1i    = findcells(RING,'FamName','O1');
o2i    = findcells(RING,'FamName','O2');
o3i    = findcells(RING,'FamName','O3');
oi       = sort([o1i o2i o3i]); ocol=[200 255 50]/255;

% ---------------
% AT bend indices
% ---------------
d1i  = findcells(RING,'FamName','D1');
d2i  = findcells(RING,'FamName','D2');
d3i  = findcells(RING,'FamName','D3');
%dmi  = findcells(RING,'FamName','Dm');
%di   = sort([d1i d2i d3i dmi]); dcol = [100 100 100]/255; % color codes for plotting
di   = sort([d1i d2i d3i]); dcol = [100 100 100]/255; % color codes for plotting

% ------------------------------------------------------
% create the dipole groupings
% dipoles per achromat (for whole bloc moves):
% DM(1)...D1(1)...D2(1)...D3(1)...D2(2)...D1(2)...DM(2)
% D*i{achromat,element}: indices of the 12 bending slices 
% ------------------------------------------------------
for iA =1:nA % NB: 20 achromats for the full ring 
    nbloc1 = numel(findcells(RING,'FamName','D1'))/12/nA;
    for i = 1:nbloc1
        D1i{iA,i} = d1i(1+(i-1)*12+(iA-1)*12*nbloc1:12+(i-1)*12+(iA-1)*12*nbloc1);
    end
    nbloc2 = numel(findcells(RING,'FamName','D2'))/12/nA;
    for i = 1:nbloc2
        D2i{iA,i} = d2i(1+(i-1)*12+(iA-1)*12*nbloc2:12+(i-1)*12+(iA-1)*12*nbloc2);
    end  
    nbloc3 = numel(findcells(RING,'FamName','D3'))/12/nA;
    for i = 1:nbloc3
        D3i{iA,i} = d3i(1+(i-1)*12+(iA-1)*12*nbloc3:12+(i-1)*12+(iA-1)*12*nbloc3);
    end    
end

% -----------------
% AT magnet indices
% -----------------
magi = unique(sort([qi ri si oi di]));


% -----------------------------------------
% create "girders" (7 girders for achromat)
% -----------------------------------------
mGs = findcells(RING,'FamName','GRDs');
mGe = findcells(RING,'FamName','GRDe');

for jA = 1:nA % 20 achromats for the full ring
    for iG = 1:7
        Gi{jA,iG} = [mGs(iG+(jA-1)*7):mGe(iG+(jA-1)*7)];
        for ii = Gi{jA,iG}
            RING{ii}.GRD=[jA, iG]; % describes which girder this element belongs to ...
        end
    end
end


% -----------------
% AT BPM indices
% -----------------
bpmi = unique(sort(findcells(RING,'FamName','BPM')));

% -----------------------
% move girders 
% -----------------------

for iA = 1:nA
    for iG = 1:7
        % -----------------------------------------
        % identify magnets within the girder limits
        % -----------------------------------------
        sel_mag  = magi(magi>Gi{iA,iG}(1)&magi<Gi{iA,iG}(end));
        % ---------------------------------------------------------------
        % identify BPMs within the girder limits >> introduce BPM offsets
        % ---------------------------------------------------------------
        sel_bpm  = bpmi(bpmi>Gi{iA,iG}(1)&bpmi<Gi{iA,iG}(end));

        % ---------------------------------------------------------------
        % identify CMs within the girder limits >> introduce CM offsets
        % ---------------------------------------------------------------
        % sel_cm  = cmi(cmi>Gi{iA,iG}(1)&cmi<Gi{iA,iG}(end));

        % -----------------------------------------
        % identify quads within the girder limits
        % -----------------------------------------
        % sel_quad = qi(qi>Gi{iA,iG}(1)&qi<Gi{iA,iG}(end));

        % -----------------------------------------
        % identify sexts within the girder limits
        % -----------------------------------------
        % sel_sext = si(si>Gi{iA,iG}(1)&si<Gi{iA,iG}(end));

        % -----------------------------------------
        % identify octs within the girder limits
        % -----------------------------------------
        % sel_oct = oi(oi>Gi{iA,iG}(1)&oi<Gi{iA,iG}(end));

        % -----------------------------------------
        % identify bends within the girder limits
        % -----------------------------------------
        % sel_bend = di(di>Gi{iA,iG}(1)&di<Gi{iA,iG}(end));
    
        % -----------------------------------------
        % move girder ...
        % -----------------------------------------
        sel_all = sort(unique([sel_mag sel_bpm]));

        RING = move_grd(RING, Gi{iA,iG}, eGr{iG});%
    end
end

% -----------------------
% move the single magnets
% -----------------------
%for iA = 1:nA
    RING = move_mag(RING, qi, eQ,'Q');%
    RING = chgrad_mag(RING,qi,eQ,'Q');
    RING = move_mag(RING, ri, eQ,'R');%
    RING = chgrad_mag(RING,ri,eQ,'R');
%    RING = move_mag(RING, di, eD);%
    RING = move_mag(RING, si, eS,'S');%
    RING = chgrad_mag(RING,si,eS,'S');
    RING = move_mag(RING, oi, eO,'O');%  
    RING = chgrad_mag(RING,oi,eO,'O');
%end

% -----------------------
% move blocs of dipoles
% -----------------------
for iA = 1:nA
    %RING = move_mag_bloc(RING, DMi{iA,1}, eD);%    
    
    for ii = 1:size(D1i,2) %numel(D1i)
        RING = move_mag_bloc(RING, D1i{iA,ii}, eD);%
        %RING = chgrad_mag_bloc(RING, D1i{iA,ii}, eD);%
    end
    for ii = 1:size(D2i,2) %numel(D2i)
        RING = move_mag_bloc(RING, D2i{iA,ii}, eD);%
        %RING = chgrad_mag_bloc(RING, D2i{iA,ii}, eD);%
    end
    for ii = 1:size(D3i,2) %numel(D3i)
        RING = move_mag_bloc(RING, D3i{iA,ii}, eD);%
        %RING = chgrad_mag_bloc(RING, D3i{iA,ii}, eD);%
    end
end

% -------------------
% plot the result ...
% -------------------
for i= 1:length(RING)
    ddx(i)  = -RING{i}.T1(1); ddy(i) = -RING{i}.T1(3);
    ddphi(i) = -asin(RING{i}.R1(1,3));
end

% -----------------------------
% magnet spatial mis-alignments 
% -----------------------------
figure(2); clf
subplot(3,1,1); hold on
plot(spos,ddx,'.'); box on; grid on
plot(spos(di),ddx(di),'o','color',dcol,'MarkerFaceColor',dcol)
plot(spos(qi),ddx(qi),'o','color',qcol,'MarkerFaceColor',qcol)
plot(spos(ri),ddx(ri),'o','color',rcol,'MarkerFaceColor',rcol)
plot(spos(si),ddx(si),'o','color',scol,'MarkerFaceColor',scol)
plot(spos(oi),ddx(oi),'o','color',ocol,'MarkerFaceColor',ocol)
plot(spos(bpmi),ddx(bpmi),'sq','color','k','MarkerFaceColor','k')
Xbpm= ddx(bpmi); Ybpm=ddy(bpmi); 

title('horizontal plane - sway/yaw'); xlabel('S (m)'); ylabel('DX (m)')
axis([0 528/5 -400e-6 400e-6])

subplot(3,1,2); hold on
plot(spos,ddy,'.'); box on; grid on
plot(spos(di),ddy(di),'o','color',dcol,'MarkerFaceColor',dcol)
plot(spos(qi),ddy(qi),'o','color',qcol,'MarkerFaceColor',qcol)
plot(spos(ri),ddy(ri),'o','color',rcol,'MarkerFaceColor',rcol)
plot(spos(si),ddy(si),'o','color',scol,'MarkerFaceColor',scol)
plot(spos(oi),ddy(oi),'o','color',ocol,'MarkerFaceColor',ocol)
plot(spos(bpmi),ddy(bpmi),'sq','color','k','MarkerFaceColor','k')
title('vertical plane - heave/pitch'); xlabel('S (m)'); ylabel('DY (m)')
axis([0 528/5 -400e-6 400e-6])

subplot(3,1,3); hold on
plot(spos,ddphi,'.'); box on; grid on
plot(spos(di),ddphi(di),'o','color',dcol,'MarkerFaceColor',dcol)
plot(spos(qi),ddphi(qi),'o','color',qcol,'MarkerFaceColor',qcol)
plot(spos(ri),ddphi(ri),'o','color',rcol,'MarkerFaceColor',rcol)
plot(spos(si),ddphi(si),'o','color',scol,'MarkerFaceColor',scol)
plot(spos(oi),ddphi(oi),'o','color',ocol,'MarkerFaceColor',ocol)
title('roll'); xlabel('S (m)'); ylabel('\Phi (rad)')
axis([0 528/5 -200e-6 200e-6])

% -----------------------------
% magnet gradient errors 
% -----------------------------
figure(3); clf
subplot(3,1,1); hold on
title('quadrupole components')
sq    = spos(qi); 
for i = 1:length(qi)
    lq(i)  = RING{qi(i)}.Length;
    kq(i)  = RING{qi(i)}.PolynomB(2); 
end
sq = sq+lq/2; 
stem(sq, kq,'linewidth',3,'marker','none')

sr    = spos(ri);
for i = 1:length(ri)
    lr(i)  = RING{ri(i)}.Length;
    kr(i)  = RING{ri(i)}.PolynomB(2); 
end
sr = sr + lr/2; 
stem(sr, kr,'linewidth',3,'marker','none')
axis([0 528/5 -7 7]); xlabel('S (m)'); ylabel('K (m-2)')

sd    = spos(di);
for i = 1:length(di)
    ld(i)  = RING{di(i)}.Length;
    th(i)  = RING{di(i)}.BendingAngle;
    kd(i)  = RING{di(i)}.PolynomB(2); 
end
sd = sd + ld/2;
stem(sd, kd,'linewidth',3,'marker','none')
Bdip = th./ld * 3/.2998; 

subplot(3,1,2); hold on
title('sextupole/octupole components')
ss    = spos(si);
for i = 1:length(si)
    ls(i)   = RING{si(i)}.Length;
    k2s(i)  = RING{si(i)}.PolynomB(3); 
end
ss = ss + ls/2; yyaxis left;
stem(ss, k2s,'linewidth',3,'marker','none');
axis([0 528/5 -300 300]); ylabel('K (m-3)')
so    = spos(oi);
for i = 1:length(oi)
    lo(i)   = RING{oi(i)}.Length;
    k3o(i)  = RING{oi(i)}.PolynomB(4); 
end
so = so + lo/2; yyaxis right;
stem(so, k3o,'linewidth',3,'marker','none')
axis([0 528/5 -3000 3000]); xlabel('S (m)'); ylabel('K (m-4)')

subplot(3,1,3); hold on
title('dipole components')

%stem(sd, Bdip,'marker','none','LineWidth',3)
for i=1:numel(sd)
    b(i) = bar(sd(i), Bdip(i), 'BarWidth',ld(i),'Facecolor','cyan');
end
axis([0 528/5 -0.1 0.7]); xlabel('S (m)'); ylabel('B (T)')

RINGs = RING; 




% -----------------
% ORBIT CORRECTIONS
% -----------------

% ensure we have 6D 
% [RING,radelemIndex,cavitiesIndex,energy] = atenable_6d(RING);
% % % % if 1 == 0
% % % % invRMx = pinv(ModelRM.OrbHCor{1});
% % % % invRMy = pinv(ModelRM.OrbVCor{3});
% % % % %orbit  = findorbit6(RING,1:length(RING)+1);
% % % % 
% % % % for io =1:4
% % % % orbit  = findorbit4(RING,1:length(RING)+1);
% % % % X      = orbit(1,:);
% % % % Y      = orbit(3,:);
% % % % for i=1:numel(bpmi)
% % % %     Xbpm(i)      = orbit(1,bpmi(i)) - (RING{bpmi(i)}.T1(1));
% % % %     Ybpm(i)      = orbit(3,bpmi(i)) - (RING{bpmi(i)}.T1(3));
% % % % end
% % % % DCMh = invRMx*Xbpm';
% % % % DCMv = invRMy*Ybpm';
% % % % 
% % % % indHCor = findcells(RING,'FamName','CMh');
% % % % indVCor = findcells(RING,'FamName','CMv');
% % % % 
% % % % for i=1:numel(indHCor)
% % % %     RING{indHCor(i)}.KickAngle(1) = RING{indHCor(i)}.KickAngle(1) -DCMh(i);
% % % % end
% % % % for i=1:numel(indVCor)
% % % %     RING{indVCor(i)}.KickAngle(2) = RING{indVCor(i)}.KickAngle(2) -DCMv(i);
% % % % end

% % % % end


%orbitc  = findorbit6(RING,1:length(RING)+1);
% % % % orbitc  = findorbit4(RING,1:length(RING)+1);
% % % % 
% % % % Xc      = orbitc(1,:);
% % % % Yc      = orbitc(3,:);
% % % % 
% % % % 
% % % % 
% % % % subplot(3,1,1); hold on
% % % % plot(spos,X(1:end-1)); 
% % % % plot(spos,Xc(1:end-1));
% % % % subplot(3,1,2); hold on
% % % % plot(spos,Y(1:end-1)); 
% % % % plot(spos,Yc(1:end-1));
% % % % end
end
% % % % % return
% % % % % 
function RING = move_grd(RING, Gi, eG)
%
% entryP, exitP 
%
% sway-yaw
entryP = Gi(1); exitP = Gi(end); 
ss     = findspos(RING,Gi);
s1     = findspos(RING,Gi(1));
s2     = findspos(RING,Gi(end));
sm     = (s1+s2)/2; % pivotal centre 

% sway/yaw
sway   = eG(1)* randn(); %sway = 123e-6; 
yaw    = eG(3)* randn(); %yaw  = 20e-6; 
dx     = sway + yaw*(ss-sm);
% heave/pitch
heave  = eG(2)* randn(); %heave = -123e-6;
pitch  = eG(4)* randn(); %pitch = -20e-6; 
dy     = heave + pitch*(ss-sm);
%RING = ataddshift(RING, Gi, dx, dy);

RING = atsetshift(RING, Gi, dx, dy,'RelativeShift');

% roll
dphi   = eG(5)* randn(); %dphi   = 30e-6; 

RING = atsettilt(RING, Gi, dphi,'RelativeTilt'); 

end

% ----------------------------------
% change gradient individual magnets
% ----------------------------------

function RING = chgrad_mag(RING, mi, me, mtyp)

for i = 1:numel(mi)
    if strcmpi(mtyp,'Q')||strcmpi(mtyp,'R')
        k0(i)  = RING{mi(i)}.PolynomB(2);
        k(i)   = k0(i) .* (1 + me(1)*trunc_randn(1,2));    
        RING{mi(i)}.PolynomB(2) = k(i);     
    elseif strcmpi(mtyp,'S')
        k0(i)  = RING{mi(i)}.PolynomB(3);           
        k(i)   = k0(i) .* (1 + me(1)*trunc_randn(1,2));
        RING{mi(i)}.PolynomB(3) = k(i);   
    elseif strcmpi(mtyp,'O')
        k0(i)  = RING{mi(i)}.PolynomB(4);                      
        k(i)   = k0(i) .* (1 + me(1)*trunc_randn(1,2)); 
        RING{mi(i)}.PolynomB(4) = k(i);   
    end
end
end

% -----------------------
% move individual magnets
% -----------------------

function RING = move_mag(RING, mi, me, mtyp)
        %
        % input: mi, magnet index / me: magnet error
        %
% %         dx     = me(2) * ones(numel(mi),1); % trunc_randn(numel(mi),1);
% %         dy     = me(3) * ones(numel(mi),1); % trunc_randn(numel(mi),1);
        dx     = me(2) * trunc_randn(numel(mi),2);
        dy     = me(3) * trunc_randn(numel(mi),2);
%        RING = ataddshift(RING, mi, dx, dy);
        RING   = atsetshift(RING, mi, dx, dy,'RelativeShift');
                
        dphi   = me(4) * trunc_randn(numel(mi),2); %trunc_randn(1,1) * ones(numel(mi),1);
        RING   = atsettilt(RING, mi, dphi,'RelativeTilt');

end


% ----------------------------------
% change gradient in bloc magnets
% ----------------------------------

function RING = chgrad_mag_bloc(RING, mi, me)

%
% define slices weight
%
% for i = 1:numel(mi)
%     w(i) = RING{mi(i)}.PolynomB(2)/RING{mi(7)}.PolynomB(2);
% end

%
% re-define quad gradients ... 
%
me_bloc = 0; % me(1)*trunc_randn(1,1);

k0 = 0.0; 
for i = 1:numel(mi)          
    k0 = RING{mi(i)}.PolynomB(2);
    RING{mi(i)}.PolynomB(2) = k0; %  * (1.0 + 0.0);     
end

%
% re-define dipole component ... 
%
if 1 == 0
for i = 1:numel(mi)          
    ba0(i)  = RING{mi(i)}.BendingAngle;
    ba(i)   = RING{mi(i)}.PolynomB(1) + ba0(i) .* me_bloc;    
    RING{mi(i)}.PolynomB(1) = ba(i);     
end
end


end
        
% -----------------------
% move magnet blocs
% -----------------------
function RING = move_mag_bloc(RING, mi, me)
        %
        % input: mi, magnet index / me: magnet error
        %
% %         dx     = me(2) * ones(numel(mi),1); % trunc_randn(numel(mi),1);
% %         dy     = me(3) * ones(numel(mi),1); % trunc_randn(numel(mi),1);
        dx     = me(2) * trunc_randn(1,2) * ones(numel(mi),1);
        dy     = me(3) * trunc_randn(1,2) * ones(numel(mi),1);
%        RING = ataddshift(RING, mi, dx, dy);
        RING = atsetshift(RING, mi, dx, dy,'RelativeShift');
        dphi   = me(4) * trunc_randn(1,2) * ones(numel(mi),1);
        RING = atsettilt(RING, mi, dphi,'RelativeTilt');

end
% % % % % 
% % % % % 
% % % % % qfcounter=0;
% % % % % qfmcounter=0;
% % % % % qfecounter=0;
% % % % % qdecounter=0;
% % % % % 
% % % % % sficounter=0;
% % % % % sfmcounter=0;
% % % % % sfocounter=0;
% % % % % sdecounter=0;
% % % % % sdcounter =0;
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % ---------------------------------------
% % % % % % no-girder move, only individual errors
% % % % % % ---------------------------------------
% % % % % if strcmpi(girder_move_type,'nomove') % if no-move, then one may want to spoil all the magnet gradients/single positions
% % % % %     sel_quad = iquad;
% % % % %     sel_sext = isext;
% % % % %     sel_oct  =  ioct;
% % % % %     if ~isempty(sel_quad)
% % % % %         [~, qfcounter, qfmcounter, qfecounter, qdecounter] = spoil_quad(eQ, sel_quad, qfcounter, qfmcounter, qfecounter, qdecounter);
% % % % %     end
% % % % %     if ~isempty(sel_sext)
% % % % %        [~, sficounter, sfmcounter, sfocounter, sdecounter, sdcounter ] = spoil_sext(eS, sel_sext, ...
% % % % %                     sficounter, sfmcounter, sfocounter, sdecounter, sdcounter); 
% % % % %     end
% % % % %     if ~isempty(sel_oct)
% % % % %         spoil_oct(eO, sel_oct);
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % % ------------------------------------------------------------------------
% % % % % % move the girders in a random or controlled mode, plus individual errors
% % % % % % ------------------------------------------------------------------------
% % % % % if ~strcmpi(girder_move_type,'nomove')
% % % % %     for iA = 1:nA   % achromats
% % % % %         for iG = 1:7   % girders per achromat
% % % % %             disp(['processing girder (' num2str(iA) ',' num2str(iG) ')'])
% % % % %             grd_num = (iA-1)*7 + iG;
% % % % %             iGRDs = findcells(RING,'FamName','GRDs');
% % % % %             iGRDe = findcells(RING,'FamName','GRDe');
% % % % % 
% % % % %             igrd_lim = [iGRDs(grd_num) iGRDe(grd_num)];
% % % % % 
% % % % %             % -----------------------------------------
% % % % %             % identify magnets within the girder limits
% % % % %             % -----------------------------------------
% % % % %             sel_mag  = magi(magi>igrd_lim(1)&magi<igrd_lim(2));
% % % % % 
% % % % %             % ---------------------------------------------------------------
% % % % %             % identify BPMs within the girder limits >> introduce BPM offsets
% % % % %             % ---------------------------------------------------------------
% % % % %             sel_bpm  = bpmi(bpmi>igrd_lim(1)&bpmi<igrd_lim(2));
% % % % % 
% % % % %             % -----------------------------------------
% % % % %             % identify quads within the girder limits
% % % % %             % -----------------------------------------
% % % % %             sel_quad = iquad(iquad>igrd_lim(1)&iquad<igrd_lim(2));
% % % % % 
% % % % %             % -----------------------------------------
% % % % %             % identify sexts within the girder limits
% % % % %             % -----------------------------------------
% % % % %             sel_sext = isext(isext>igrd_lim(1)&isext<igrd_lim(2));
% % % % % 
% % % % %             % -----------------------------------------
% % % % %             % identify octs within the girder limits
% % % % %             % -----------------------------------------
% % % % %             sel_oct = ioct(ioct>igrd_lim(1)&ioct<igrd_lim(2));
% % % % % 
% % % % %             % -----------------------------------------
% % % % %             % identify bends within the girder limits
% % % % %             % -----------------------------------------
% % % % %             sel_bend = ibend(ibend>igrd_lim(1)&ibend<igrd_lim(2));
% % % % % 
% % % % %             % ------------------------------------------------------
% % % % %             % spoil the gradients and introduce individual magnet
% % % % %             % displacements of the families within a selected girder
% % % % %             % ------------------------------------------------------
% % % % % 
% % % % %             if ~isempty(sel_quad)
% % % % %                 [~, qfcounter, qfmcounter, qfecounter, qdecounter ] = spoil_quad(eQ, sel_quad, ...
% % % % %                     qfcounter, qfmcounter, qfecounter, qdecounter);
% % % % %             end
% % % % %             if ~isempty(sel_sext)
% % % % %                 [~, sficounter, sfmcounter, sfocounter, sdecounter, sdcounter ] = spoil_sext(eS, sel_sext, ...
% % % % %                     sficounter, sfmcounter, sfocounter, sdecounter, sdcounter);
% % % % % %                spoil_sext(eS, sel_sext);
% % % % %             end
% % % % %             if ~isempty(sel_oct)
% % % % %                 spoil_oct(eO, sel_oct);
% % % % %             end
% % % % % 
% % % % %             % ---------------
% % % % %             % move the girder
% % % % %             % ---------------
% % % % %             move_grd(sel_mag, sel_bpm, bpmi, igrd_lim, grd_num, gr, gt, girder_move_type);
% % % % % 
% % % % %         end
% % % % %     end
% % % % % end
% % % % % 
% % % % % %
% % % % % % ad HOC test
% % % % % %
% % % % % % setsp('QF',  4.038136152,'physics','simulator')
% % % % % % setsp('QFM', 3.781542995,'physics','simulator')
% % % % % % setsp('QFE', 3.77397990402,'physics','simulator')
% % % % % % setsp('QDE',-2.4986556740,'physics','simulator')
% % % % % 
% % % % % RING = RING; 
% % % % % 
% % % % % %
% % % % % % post checks
% % % % % %
% % % % % graphic = 'true';
% % % % % oc      = OrbitCorrectionFlag; 
% % % % % fakeBBA = 'true';
% % % % % 
% % % % % 
% % % % % %%% 430e_VC_GRD ... RM = load('/gpfs/offline1/staff/common/marapo/atcollabMML/machine/MAXIV/R3OpsData/ProductionOptics-2022-RP1/GoldenBPMResp_R3_ProductionOptics2022_RP1.mat');
% % % % % %%% RM = load('/gpfs/offline1/staff/common/marapo/atcollabMML/machine/MAXIV/R3OpsData/ProductionOptics-2020-RP1/GoldenBPMResp_R3_ProductionOptics2020_RP1.mat')
% % % % % RM = load('/gpfs/offline1/staff/common/marapo/atcollabMML/machine/MAXIV/R3OpsData/ProductionOptics-2021-RP1/GoldenBPMResp_R3_ProductionOptics2021_RP1.mat');
% % % % % BPMdevlist = RM.Rmat(1,1).Monitor.DeviceList;
% % % % % sbpm       = getspos('BPMx',BPMdevlist);
% % % % % BBA_errorx = 0e-6; BBA_errory  = 0e-6; 
% % % % % GoalOrbitX = getoffset('BPMx', BPMdevlist); 
% % % % % GoalOrbitY = getoffset('BPMy', BPMdevlist);
% % % % % if fakeBBA
% % % % %     GoalOrbitX = GoalOrbitX + randn_t(length(GoalOrbitX),1,2) * BBA_errorx;
% % % % %     GoalOrbitY = GoalOrbitY + randn_t(length(GoalOrbitY),1,2) * BBA_errory;
% % % % % end     
% % % % % 
% % % % % bpmoffx = GoalOrbitX;
% % % % % bpmoffy = GoalOrbitY;
% % % % % 
% % % % % 
% % % % % if OrbitCorrectionFlag
% % % % % %for k = 1:1
% % % % % ErrorFlagx=1; ErrorFlagy=1; neigeny = 110; neigenx = 160; loopn=0;
% % % % % while ErrorFlagx==1 || ErrorFlagy ==1
% % % % %     loopn = loopn+1;
% % % % %     disp(['loop n = ' num2str(loopn)])
% % % % %     [OCS, OCS0, V, S, ErrorFlagx]  = setorbitMA({GoalOrbitX}, {getx('physics','struct',BPMdevlist)},{getsp('HCM','physics','struct')},2,neigenx,'CorrectorGain', 0.85);
% % % % %     if ErrorFlagx == 1
% % % % %         setsp('HCM',0,'simulator','physics');
% % % % %         neigenx = neigenx - 10;
% % % % %     end
% % % % %     [OCS, OCS0, V, S, ErrorFlagy]  = setorbitMA({GoalOrbitY}, {gety('physics','struct',BPMdevlist)},{getsp('VCM','physics','struct')},2,neigeny,'CorrectorGain', 0.85);
% % % % %     if ErrorFlagy == 1
% % % % %         setsp('VCM',0,'simulator','physics');
% % % % %         neigeny = neigeny - 10;
% % % % %     end
% % % % % 
% % % % %     if neigenx <=0 || neigeny <= 0
% % % % %         ErrorFlagx = 0;
% % % % %         ErrorFlagy = 0;
% % % % %     end
% % % % % 
% % % % % end
% % % % % 
% % % % % ErrorFlagx=1; ErrorFlagy=1; neigeny = 140; neigenx = 180; loopn=0;
% % % % % while ErrorFlagx==1 || ErrorFlagy ==1
% % % % %     loopn = loopn+1;
% % % % %     disp(['loop n = ' num2str(loopn)])
% % % % %     [OCS, OCS0, V, S, ErrorFlagx]  = setorbitMA({GoalOrbitX}, {getx('physics','struct',BPMdevlist)},{getsp('HCM','physics','struct')},2,neigenx,'CorrectorGain', 0.90);
% % % % %     if ErrorFlagx == 1
% % % % %         setsp('HCM',0,'simulator','physics');
% % % % %         neigenx = neigenx - 10;
% % % % %     end
% % % % %     [OCS, OCS0, V, S, ErrorFlagy]  = setorbitMA({GoalOrbitY}, {gety('physics','struct',BPMdevlist)},{getsp('VCM','physics','struct')},2,neigeny,'CorrectorGain', 0.90);
% % % % %     if ErrorFlagy == 1
% % % % %         setsp('VCM',0,'simulator','physics');
% % % % %         neigeny = neigeny - 10;
% % % % %     end
% % % % %     if neigenx <=0 || neigeny <= 0
% % % % %         ErrorFlagx = 0;
% % % % %         ErrorFlagy = 0;
% % % % %     end
% % % % % 
% % % % % end
% % % % % 
% % % % % ErrorFlagx=1; ErrorFlagy=1; neigeny = 183; neigenx = 200; loopn=0;
% % % % % while ErrorFlagx==1 || ErrorFlagy ==1
% % % % %     loopn = loopn+1;
% % % % %     disp(['loop n = ' num2str(loopn)])
% % % % %     [OCS, OCS0, V, S, ErrorFlagx]  = setorbitMA({GoalOrbitX}, {getx('physics','struct',BPMdevlist)},{getsp('HCM','physics','struct')},3,neigenx,'CorrectorGain', 0.95);
% % % % %     if ErrorFlagx == 1
% % % % %         setsp('HCM',0,'simulator','physics');
% % % % %         neigenx = neigenx - 10;
% % % % %     end
% % % % %     [OCS, OCS0, V, S, ErrorFlagy]  = setorbitMA({GoalOrbitY}, {gety('physics','struct',BPMdevlist)},{getsp('VCM','physics','struct')},3,neigeny,'CorrectorGain', 0.95);
% % % % %     if ErrorFlagy == 1
% % % % %         setsp('VCM',0,'simulator','physics');
% % % % %         neigeny = neigeny - 10;
% % % % %     end
% % % % %     if neigenx <=0 || neigeny <= 0
% % % % %         ErrorFlagx = 0;
% % % % %         ErrorFlagy = 0;
% % % % %         return
% % % % %     end
% % % % % 
% % % % % end
% % % % % end
% % % % % 
% % % % % % % % % % for k = 1:1
% % % % % % % % % %     [OCS, OCS0, V, S, ErrorFlag]  = setorbitMA({GoalOrbitY}, {gety('physics','struct',BPMdevlist)},{getsp('VCM','physics','struct')},3,130,'CorrectorGain', 0.90);
% % % % % % % % % %     [OCS, OCS0, V, S, ErrorFlag]  = setorbitMA({GoalOrbitX}, {getx('physics','struct',BPMdevlist)},{getsp('HCM','physics','struct')},3,155,'CorrectorGain', 0.90);
% % % % % % % % % % end
% % % % % % % % % % for k = 1:1
% % % % % % % % % %     [OCS, OCS0, V, S, ErrorFlag]  = setorbitMA({GoalOrbitY}, {gety('physics','struct',BPMdevlist)},{getsp('VCM','physics','struct')},3,183,'CorrectorGain', 0.95);
% % % % % % % % % %     [OCS, OCS0, V, S, ErrorFlag]  = setorbitMA({GoalOrbitX}, {getx('physics','struct',BPMdevlist)},{getsp('HCM','physics','struct')},3,200,'CorrectorGain', 0.95);
% % % % % % % % % % end
% % % % % [oc s] = getorbit_everywhere_maxiv;
% % % % % HCM = getsp('HCM','physics','simulator'); shcm = getspos('HCM');
% % % % % VCM = getsp('VCM','physics','simulator'); svcm = getspos('VCM');
% % % % % 
% % % % % RING = RING; 
% % % % % 
% % % % % % --------------------------------------------
% % % % % % plot the girder mis.alignments / CM / orbits
% % % % % % --------------------------------------------
% % % % % 
% % % % % if graphic
% % % % % 
% % % % %     figure(1); clf; hold on; tiledlayout(3,1)
% % % % %     ax1 = nexttile; 
% % % % %     ax2 = nexttile; 
% % % % %     ax3 = nexttile; 
% % % % %     linkaxes([ax1 ax2 ax3],'x')
% % % % %     smag      = findspos(RING,magi); % magnet positions (entry point)
% % % % %     for k = 1:length(magi)
% % % % %         lmag(k)    = RING{magi(k)}.Length;
% % % % %         smag(k)    = smag(k) + lmag(k)/2; % magnet positions (magnet centre)
% % % % %     end
% % % % % 
% % % % % 
% % % % %     for k = 1:length(smag)
% % % % %         dx_mag(k) = RING{magi(k)}.T2(1); dy_mag(k) = RING{magi(k)}.T2(3);
% % % % %     end
% % % % %     
% % % % %     %subplot(3,1,2); hold on
% % % % %     axes(ax2)
% % % % %     plot(smag, dx_mag,'o'); hold on; plot(smag,dy_mag,'o'); 
% % % % %     axis([0 528 -6e-4 6e-4]); box on; grid on
% % % % %     legend('H-plane','V-plane'); title('girder positions')
% % % % %     ylabel('orbit (m)')
% % % % %     
% % % % %     
% % % % %     figure(1); hold on; axes(ax2)
% % % % % plot(sbpm, bpmoffx,'sb','MarkerFaceColor','b'); hold on
% % % % % plot(sbpm, bpmoffy,'sr','MarkerFaceColor','r');
% % % % % legend('H-plane','V-plane','BPMx','BPMy')
% % % % % 
% % % % % 
% % % % % figure(1); hold on;  axes(ax2)
% % % % % %subplot(3,1,2); hold on 
% % % % % plot(s, oc(1,:),'b'); plot(s, oc(3,:),'r')
% % % % % box on; grid on
% % % % % legend('H-plane','V-plane','BPMx','BPMy','cOx','cOy');
% % % % % 
% % % % % figure(1); hold on ;  axes(ax1)
% % % % % %subplot(3,1,1); hold on
% % % % % stem(shcm,HCM,'b','marker','none','linewidth',3); hold on
% % % % % stem(svcm,VCM,'r','marker','none','linewidth',3);
% % % % % box on; grid on; legend('HCM','VCM'); 
% % % % % title('corrector strength'); ylabel('CM (rad)')
% % % % % legend('HCM','VCM'); axis([0 528 -400e-6 400e-6])
% % % % % 
% % % % % 
% % % % % figure(1); hold on ;  axes(ax3)
% % % % % %subplot(3,1,3); hold on
% % % % % plot(smag, 1e6*(dx_mag-oc(1,magi)),'.b'); hold on
% % % % % plot(smag, 1e6*(dy_mag-oc(3,magi)),'.r')
% % % % % legend('resH-plane','resV-plane'); box on; grid on
% % % % % title('residuals at magnets'); ylabel('\deltaXY  (\mum)')
% % % % % axis([0 528 -350 350])
% % % % % end
% % % % % 
% % % % % 
% % % % % end
% % % % % 
% % % % % 
% % % % % function R = move_grd(sel_mag, sel_bpm, ibpm, igrd_lim, grd_num, gr, gt, misal_type)
% % % % % global RING
% % % % %     bpm_dev = family2dev('BPM',0);
% % % % %    
% % % % %     nn = mod(grd_num,7);
% % % % %     nn(nn==0) = 7; 
% % % % % 
% % % % %     switch misal_type
% % % % %         case 'random'
% % % % %             Ssway  = gr{nn}(1);
% % % % %             Sheave = gr{nn}(2);
% % % % %             Syaw   = gr{nn}(3);
% % % % %             Spitch = gr{nn}(4);
% % % % %             Sroll  = gr{nn}(4);
% % % % %             sway   = Ssway  * randn_t(1,1,2);
% % % % %             heave  = Sheave * randn_t(1,1,2);
% % % % %             yaw    = Syaw   * randn_t(1,1,2);
% % % % %             pitch  = Spitch * randn_t(1,1,2);
% % % % %             roll   = Sroll  * randn_t(1,1,2);
% % % % %         case 'fixed'
% % % % %             Ssway  = gt{nn}(1);
% % % % %             Sheave = gt{nn}(2);
% % % % %             Syaw   = gt{nn}(3);
% % % % %             Spitch = gt{nn}(4);
% % % % %             Sroll  = gt{nn}(4);
% % % % %             sway   = Ssway  ;
% % % % %             heave  = Sheave ;
% % % % %             yaw    = Syaw   ;
% % % % %             pitch  = Spitch ;
% % % % %             roll   = Sroll  ;            
% % % % %     end
% % % % %     
% % % % %     %
% % % % %     % implement sway and heave
% % % % %     %
% % % % %     [dx0, dy0] = getshift(sel_mag);
% % % % %     setshift(sel_mag, dx0 + sway, dy0 + heave);
% % % % % 
% % % % %     [~, b] = intersect(ibpm, sel_bpm);
% % % % %     xoff0  = getoffset('BPMx',bpm_dev(b,:));
% % % % %     yoff0  = getoffset('BPMy',bpm_dev(b,:));
% % % % %     setoffset('BPMx',xoff0+sway,bpm_dev(b,:));  % -sway  because if you move by Dx you will see the orbit at -Dx 
% % % % %     setoffset('BPMy',yoff0+heave,bpm_dev(b,:)); % -heave because if you move by Dy you will see the orbit at -Dy 
% % % % %     
% % % % %     
% % % % %     %
% % % % %     % implement yaw and pitch
% % % % %     %
% % % % % 
% % % % %     % for magnets
% % % % %     sgrd_lim = findspos(RING, igrd_lim);
% % % % %     smid     = mean(sgrd_lim); % mid point of the girder (to be reviewed)
% % % % %     s = findspos(RING, sel_mag);
% % % % %     dx_yaw   = yaw * (s - smid); 
% % % % %     dy_yaw   = zeros(size(dx_yaw));
% % % % %     dy_pitch = pitch * (s - smid); 
% % % % %     dx_pitch = zeros(size(dy_pitch));
% % % % %    
% % % % %     [dx0, dy0] = getshift(sel_mag);
% % % % %     setshift(sel_mag, dx0 + dx_yaw, dy0 + dy_pitch);
% % % % % 
% % % % %     % for BPMs 
% % % % %     s = findspos(RING, sel_bpm);
% % % % %     dx_yaw   = yaw * (s - smid); 
% % % % %     dy_yaw   = zeros(size(dx_yaw));
% % % % %     dy_pitch = pitch * (s - smid); 
% % % % %     dx_pitch = zeros(size(dy_pitch));
% % % % % 
% % % % %     [~, b] = intersect(ibpm, sel_bpm);
% % % % %     xoff0  = getoffset('BPMx',bpm_dev(b,:));
% % % % %     yoff0  = getoffset('BPMy',bpm_dev(b,:));
% % % % %     setoffset('BPMx',xoff0+dx_yaw',bpm_dev(b,:));  % -sway  because if you move by Dx you will see the orbit at -Dx 
% % % % %     setoffset('BPMy',yoff0+dy_pitch',bpm_dev(b,:)); % -heave because if you move by Dy you will see the orbit at -Dy 
% % % % % 
% % % % %     bpmoffx = getoffset('BPMx',bpm_dev(b,:));
% % % % %     bpmoffy = getoffset('BPMy',bpm_dev(b,:));
% % % % % end
% % % % % 
% % % % % function [R , qfcounter, qfmcounter, qfecounter, qdecounter]= spoil_quad(eQ, sel_quad, qfcounter, qfmcounter, qfecounter, qdecounter)
% % % % % global RING
% % % % % 
% % % % % % ------------------------------------------------------------------
% % % % % % introduce  gradient errors - use AT to introduce individual errors
% % % % % % ------------------------------------------------------------------
% % % % % 
% % % % % 
% % % % % %
% % % % % % load external errors
% % % % % %
% % % % % exterr = load('QUADERRORS_elegant_s1004.mat');
% % % % % %exterr = load('SEXTERRORS_elegant_s1005.mat');
% % % % % 
% % % % % if 1 == 1 % exclude quads
% % % % % exterr = load('QUADERRORS_elegant_s1004.mat');
% % % % % 
% % % % % for i=1:length(sel_quad)
% % % % %     k0(i) = RING{sel_quad(i)}.PolynomB(2);
% % % % %     if strcmpi(RING{sel_quad(i)}.FamName,'qf')
% % % % %         qfcounter=qfcounter+1;
% % % % %         k(i)  = k0(i) .* (1 + exterr.QF.FSE(qfcounter));
% % % % %     elseif strcmpi(RING{sel_quad(i)}.FamName,'qfm')
% % % % %         qfmcounter=qfmcounter+1;
% % % % %         k(i)  = k0(i) .* (1 + exterr.QFM.FSE(qfmcounter));
% % % % %     elseif strcmpi(RING{sel_quad(i)}.FamName,'qfend')
% % % % %         qfecounter=qfecounter+1;
% % % % %         k(i)  = k0(i) .* (1 + exterr.QFEND.FSE(qfecounter));
% % % % %     elseif strcmpi(RING{sel_quad(i)}.FamName,'qdend')
% % % % %         qdecounter=qdecounter+1;
% % % % %         k(i)  = k0(i) .* (1 + exterr.QDEND.FSE(qdecounter));
% % % % %     end
% % % % %     RING{sel_quad(i)}.PolynomB(2) = k(i);
% % % % %     RING{sel_quad(i)}.K           = k(i);
% % % % % end
% % % % % end
% % % % % 
% % % % % % -------------------------------
% % % % % % introduce magnet mis-alignments
% % % % % % -------------------------------
% % % % % dx = eQ(2).*randn_t(length(sel_quad),1,2);
% % % % % dy = eQ(3).*randn_t(length(sel_quad),1,2);
% % % % % [dx0, dy0] = getshift(sel_quad);
% % % % % setshift(sel_quad, dx0 + dx, dy0 + dy)
% % % % % 
% % % % % R = RING; 
% % % % % end
% % % % % 
% % % % % 
% % % % % function   [R, sficounter, sfmcounter, sfocounter, sdecounter, sdcounter ] = spoil_sext(eS, sel_sext, ...
% % % % %                     sficounter, sfmcounter, sfocounter, sdecounter, sdcounter);
% % % % % 
% % % % % global RING
% % % % % 
% % % % % % ------------------------------------------------------------------
% % % % % % introduce  gradient errors - use AT to introduce individual errors
% % % % % % ------------------------------------------------------------------
% % % % % % note sext are always split in two here
% % % % % if 1 == 0; % exclude sexts
% % % % % exterr = load('SEXTERRORS_elegant_s1005.mat');
% % % % % 
% % % % % for i=1:length(sel_sext)
% % % % %     k0(i) = RING{sel_sext(i)}.PolynomB(3);
% % % % % 
% % % % %     if strcmpi(RING{sel_sext(i)}.FamName,'sfi')
% % % % %         sficounter=sficounter+1;
% % % % %         sfic    = floor(0.5*sficounter) + mod(sficounter,2);
% % % % %         k(i)  = k0(i) .* (1 + exterr.SFI.FSE(sfic));
% % % % %     elseif strcmpi(RING{sel_sext(i)}.FamName,'sfm')            
% % % % %         sfmcounter=sfmcounter+1;
% % % % %         sfmc = floor(0.5*sfmcounter) + mod(sfmcounter,2);
% % % % %         k(i)  = k0(i) .* (1 + exterr.SFM.FSE(sfmc));
% % % % %     elseif strcmpi(RING{sel_sext(i)}.FamName,'sfo')
% % % % %         sfocounter= sfocounter + 1;
% % % % %         sfoc      = floor(0.5*sfocounter) + mod(sfocounter,2);
% % % % %         k(i)  = k0(i) .* (1 + exterr.SFO.FSE(sfoc));
% % % % %     elseif strcmpi(RING{sel_sext(i)}.FamName,'sdend')
% % % % %         sdecounter= sdecounter + 1;
% % % % %         sdec      = floor(0.5*sdecounter) + mod(sdecounter,2);
% % % % %         k(i)  = k0(i) .* (1 + exterr.SDE.FSE(sdec));
% % % % %     elseif strcmpi(RING{sel_sext(i)}.FamName,'sd')
% % % % %         sdcounter=sdcounter+1;
% % % % %         sdc      = floor(0.5*sdcounter) + mod(sdcounter,2);
% % % % %         k(i)  = k0(i) .* (1 + exterr.SD.FSE(sdc));
% % % % %     end
% % % % %     RING{sel_sext(i)}.PolynomB(3) = k(i);
% % % % % end
% % % % % end
% % % % % 
% % % % % % -------------------------------
% % % % % % introduce magnet mis-alignments
% % % % % % -------------------------------
% % % % % dx = eS(2).*randn_t(length(sel_sext),1,2);
% % % % % dy = eS(3).*randn_t(length(sel_sext),1,2);
% % % % % [dx0, dy0] = getshift(sel_sext);
% % % % % setshift(sel_sext, dx0 + dx, dy0 + dy)
% % % % % 
% % % % % 
% % % % % R = RING; 
% % % % % end
% % % % % 
% % % % % 
% % % % % function R = spoil_oct(eO, sel_oct)
% % % % % global RING
% % % % % 
% % % % % % ------------------------------------------------------------------
% % % % % % introduce  gradient errors - use AT to introduce individual errors
% % % % % % ------------------------------------------------------------------
% % % % % for i=1:length(sel_oct)
% % % % %     k0(i) = RING{sel_oct(i)}.PolynomB(4);
% % % % %     k(i)  = k0(i) .* (1 + eO(1).*randn_t(length(k0(i)),1,2));
% % % % %     RING{sel_oct(i)}.PolynomB(4) = k(i);
% % % % % end
% % % % % 
% % % % % % -------------------------------
% % % % % % introduce magnet mis-alignments
% % % % % % -------------------------------
% % % % % dx = eO(2).*randn_t(length(sel_oct),1,2);
% % % % % dy = eO(3).*randn_t(length(sel_oct),1,2);
% % % % % [dx0, dy0] = getshift(sel_oct);
% % % % % setshift(sel_oct, dx0 + dx, dy0 + dy)
% % % % % 
% % % % % 
% % % % % R = RING; 
% % % % % end
% % % % % 
% % % % % 
% % % % % function v = randn_t(a,b, trunc)
% % % % % outlier = ones(a,1);
% % % % % v = randn(a,b);
% % % % % while sum(outlier)>0;
% % % % %     v(find(outlier>0)) = randn(length(find(outlier>0)),b);
% % % % %     outlier = abs(v)>trunc; % trunc-sigma truncation
% % % % % end
% % % % % end