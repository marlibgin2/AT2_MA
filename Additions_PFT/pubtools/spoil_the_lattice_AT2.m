function varargout = spoil_the_lattice_AT2(varargin)%, OrbitCorrectionFlag)
% SPOIL_THE_LATTICE deployes misalignments for sliced magnets and 
% [RINGWE, RING0, MAGe, GIRe, Xbpm, Ybpm] = spoil_the_lattice_AT2(RING)
%
% INPUT
% 1. RING   -   {cell array of structs} AT2 lattice, with or without
%               errors. Note that in order to correctly assign correlated
%               errors there are some pre-requisites. In particular:
%                   a) Correlation of errors over a girder require that the
%                      girder start and end have marker elements with FamName
%                      'GS' and 'GE' in the lattice.
%                   b) Correlation of errors for sliced magnets are
%                      possible via the use of a 'MagNum' field. Elements with
%                      the same 'MagNum' are assumed to belong to the same
%                      magnet. Note that if the field does not exist for
%                      any element in the input lattice a best guess will
%                      be made. See markSlicedMagnets and getmagnetslices
%                      for details.
% 2. Flags:
%       'Display'       Applied errors will be plotted
%       'NoDisplay'     {default} No plot is generated. Note that
%                       atplot(RING,@plotMisalignments) may be used to view
%                       the errors at a later stage
%
%       'Relative'      {default} Errors are added to any already existing in the
%                       input lattice. NB! Always the case for field errors! 
%       'Absolute'      Pre-existing misalignment errors are overwritten.
%                       Field errors are not affected.
%
% OUTPUT
% 1. RINGWE -   {cell array of structs} AT2 lattice output, with new errors
%               applied.
% 2. RING0  -   {cell array of structs} AT2 lattice output, without errors.
% 
% See also markSlicedMagnets, getmagnetslices, getMagGroupsFromGirderIndex

%T = load('ModelRM.mat'); ModelRM = T.ModelRM; clear T

%% Default settings
DisplayFlag = false; % plot final results
ExpandFlag = true;
RelativeFlag = true;

%% Input handling
if nargin < 2
    nA = 1; % number of achromats
end

for n = nargin:-1:1
    if iscell(varargin{n})
        if isstruct(varargin{n}{1})
            RING = varargin{n};
            varargin(n) = [];
        end
    elseif ischar(varargin{n}) || isStringScalar(varargin{n})
        switch lower(varargin{n})
            case 'display'
                DisplayFlag = true;
                varargin(n) = [];
            case 'nodisplay'
                DisplayFlag = false;
                varargin(n) = [];
            case 'expand'
                ExpandFlag = varargin{n+1};
                varargin(n+1) = [];
                varargin(n) = [];
            case 'relative'
                RelativeFlag = true;
                varargin(n) = [];
            case 'absolute'
                RelativeFlag = false;
                varargin(n) = [];
        end
    elseif isstruct(varargin{n})
    
    end
end


% -------------------------
% Clean the input lattice
% -------------------------
% Assign zero shifts to all elements, creating the T1, T2 fields.
% This is done explicity in all dimensions as atsetshift and atsettilt
% only manipulate the transverse errors and assumptions regarding the
% input are to be avoided on general principle.
RING0 = RING;
for n = 1:numel(RING0)
    if isfield(RING0{n},'T1'), RING0{n}.T1 = zeros(6,1); RING0{n}.T2 = zeros(6,1); end
    if isfield(RING0{n},'R1'), RING0{n}.R1 = eye(6,6); RING0{n}.R1 = eye(6,6); end
end

% If misalignment errors are to be applied on a clean lattice, use the
% clean lattice instead.
if RelativeFlag == false
    RING = RING0; 
end

%% Calculate static information for quick look-up
% Get the element s-positions, which is used more or less everywhere.
allelemi = (1:length(RING))';
spos     = findspos(RING,allelemi)';

% If sliced magnets have not been marked in the lattice, do so now
% NB! markSlicedMagnets relies on all slices having the same FamName.
if isempty(findcells(RING,'MagNum'))
    RING = markSlicedMagnets(RING);
end

% Get indices for sliced magnets
% NB! Relies on MagNum field being present (set above by markSlicedMagnets)
index_magnetslices = getmagnetslices(RING);

% Get indices for elements sharing the same girder
% NB! Relies on marker elements with FamName GS and GE being present in the
% lattice to indicate girders.
index_girderelements = getMagGroupsFromGirderIndex(RING);


%% DEFINE ERRORS
% Should be given as an input to the function instead. Struct definition?

% ----------------------------------------
% single magnet error table (RMS) --> MAGe
% ----------------------------------------
%      grad(frac)   dx(um)    dy(um)
gradZero  = 0.0;  % 1 turn off/on the gradient errors
shiftZero = 0.25; % 1 turn off/on the displacement errors
rollZero  = 0;  % 1 turn off/on the roll errors
eQ = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 100e-6*rollZero ]; % quadrupole
eR = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 100e-6*rollZero ]; % reverse-bends
eS = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % sextupole
eO = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % octupole
eD = [ 5e-4*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % dipole
MAGe.eQ = eQ; MAGe.eR = eR; MAGe.eS = eS; MAGe.eO = eO; MAGe.eD = eD;

% errors.Quadrupole =
% errors.ReverseBend =
% errors.Sextupole =
% errors.Octupole =
% errors.Bend =
% errors.

% -------------------------------
% girder random error table (RMS)
% -------------------------------
% <MSj> Given the girder shape, i.e. it's longer than it's wide, the
% expectation is that roll will be harder to correctly determine. I
% therefore swapped the yaw/pitch and roll values.
%
%        sway(um) heave(um) yaw(urad) pitch(urad) roll(urad)
grdZero = 0;
eGr{1}  = [20      20       10        10          25] * 1e-6  *grdZero;
eGr{2}  = [20      20       10        10          25] * 1e-6  *grdZero;
eGr{3}  = [20      20       10        10          25] * 1e-6  *grdZero;
eGr{4}  = [20      20       10        10          25] * 1e-6  *grdZero;
eGr{5}  = [20      20       10        10          25] * 1e-6  *grdZero;
eGr{6}  = [20      20       10        10          25] * 1e-6  *grdZero;
eGr{7}  = [20      20       10        10          25] * 1e-6  *grdZero;

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


%% Deploy girder errors
for n = 1:numel(index_girderelements)

    % Depending on the girder type there may be different error magnitudes
    girderTypeIndex = mod(n-1,7)+1;
    move_grd(index_girderelements{n},eGr{girderTypeIndex});
end


%% Deploy magnet errors

for n = 1:numel(index_magnetslices)
    % Get class
    atclass = atguessclass(RING{index_magnetslices{n}(1)});

    % Get the correct error amplitudes
    switch lower(atclass)
        case 'bend'
            % NB! Errors may be different for reverse bends. Current
            % assumption is that they are similar to the quadrupole error
            % level.
            if RING{index_magnetslices{n}(1)}.BendingAngle > 0
                misalignmentError = eD;
                fieldError = eD;
            else
                misalignmentError = eQ;
                fieldError = eQ;
            end
        case 'quadrupole'
            misalignmentError = eQ;
            fieldError = eQ;
        case 'sextupole'
            misalignmentError = eS;
            fieldError = eS;
        case {'multipole','octupole'}
            % NB! atguessclass do not identify octupoles as such, only as
            % multipoles.
            misalignmentError = eO;
            fieldError = eO;
        otherwise
            fprintf('Unknown magnet class [%d - %d: %s], no error deployed.\n',index_magnetslices{n}(1),index_magnetslices{n}(end),atclass);
    end

    % Deploy errors via nested functions to speed up the execution (no need
    % to pass the RING variable to the functions)
    move_mag(index_magnetslices{n},misalignmentError);
    chgrad_mag(index_magnetslices{n},fieldError);

end

% ----------------------
% spoiling functions ...
% ----------------------
% <MSj> In the interest of speed I reworked the below into nested
% functions, to avoid the overhead of copying the lattice for each girder
% or magnet that should be moved.
% Further work is likely needed here, as atsetshift and atsettilt function
% calls also leads to a significant amount of copying.

    function move_grd(Gi, eG)
        % entryP = Gi(1); exitP = Gi(end);

        % -------------------------------------
        % define coordinates of girder elements
        % and pivotal point for yaw/pitch moves
        % -------------------------------------
        %         ss     = findspos(RING,Gi);
        s1     = spos(Gi(1));
        s2     = spos(Gi(end));
        sm     = (s1+s2)/2; % pivotal centre

        % --------------------------
        % horizontal plane: sway/yaw
        % --------------------------
        sway   = eG(1)* randn();
        yaw    = eG(3)* randn();
        dx     = sway + yaw*(spos(Gi)-sm);

        % ---------------------------
        % vertical plane: heave/pitch
        % ---------------------------
        heave  = eG(2)* randn();
        pitch  = eG(4)* randn();
        dy     = heave + pitch*(spos(Gi)-sm);

        % -------------------------------
        % move in (x,y), apply roll (phi)
        % -------------------------------
        dphi   = eG(5)* randn();

        %         RING = atsetshift(RING, Gi, dx, dy,'RelativeShift');
        %         RING = atsettilt(RING, Gi, dphi,'RelativeTilt');

        % For speed, call atshiftelem and attiltelem functions directly
        % rather than calling atsetshift and atsettilt. This avoids the
        % quite large overhead of copying the entire lattice a large number
        % of times.
        for ii = 1:length(Gi)
            RING{Gi(ii)}=atshiftelem(RING{Gi(ii)},dx(ii),dy(ii),'RelativeShift');
            RING{Gi(ii)}=attiltelem(RING{Gi(ii)},dphi,'RelativeTilt');
        end


    end

% ----------------------------------
% change gradient individual magnets
% ----------------------------------

    function chgrad_mag(mi, me)
        % -------------------------------------------
        % Affects all field components; the assumption is that they are
        % defined by the iron.
        % -------------------------------------------
  
        scalingError = me(1)*trunc_randn(1,2);

        for ii = 1:numel(mi)
            RING{mi(ii)}.PolynomB = (1 + scalingError)*RING{mi(ii)}.PolynomB;
            RING{mi(ii)}.PolynomA = (1 + scalingError)*RING{mi(ii)}.PolynomA;
         end
    end

% -----------------------
% move individual magnets
% -----------------------

    function move_mag(mi, me)
        % ------------------------------------------
        % input: mi, magnet index / me: magnet error
        % ------------------------------------------
        % Generate shifts and rolls for this single magnet
        dx     = me(2) * trunc_randn(1,2);
        dy     = me(3) * trunc_randn(1,2);
        dphi   = me(4) * trunc_randn(1,2);

        %         RING   = atsetshift(RING, mi, dx, dy,'RelativeShift');
        %         RING   = atsettilt(RING, mi, dphi,'RelativeTilt');

        % For speed, call atshiftelem and attiltelem functions directly
        % rather than calling atsetshift and atsettilt. This avoids the
        % quite large overhead of copying the entire lattice a large number
        % of times.
        for ii = 1:length(mi)
            RING{mi(ii)}=atshiftelem(RING{mi(ii)},dx,dy,'RelativeShift');
            RING{mi(ii)}=attiltelem(RING{mi(ii)},dphi,'RelativeTilt');
        end

    end




%% Assign outputs

if nargout > 0
    varargout{1} = RING;
end

if nargout > 1
    varargout{2} = RING0;
end

if nargout > 2
    varargout{3} = MAGe;
end

if nargout > 3
    varargout{4} = GIRe;
end

if nargout > 4
    % Extract errors from the lattice
    bpmi = atgetcells(RING,'FamName','BPM');
    Xbpm = -atgetfieldvalues(RING(bpmi),'T1',{1});
    Ybpm = -atgetfieldvalues(RING(bpmi),'T1',{3});
    
    varargout{5} = Xbpm;
end

if nargout > 5
    varargout{6} = Ybpm;
end


% -------------------
% plot the result ...
% -------------------

if DisplayFlag == 1
    % PLOT results

    qcol=[255 100 0]/255;
    rcol=[255 0 200]/255;
    scol=[0 255 200]/255;
    ocol=[200 255 50]/255;
    dcol = [100 100 100]/255; % color codes for plotting

    % -----------------------------
    % magnet spatial mis-alignments
    % -----------------------------
    % Extract errors from the lattice
%     for i= 1:length(RING)
%         ddx(i)  = -RING{i}.T1(1); ddy(i) = -RING{i}.T1(3);
%         ddphi(i) = -asin(RING{i}.R1(1,3)); 
%     end

    ddx     = -atgetfieldvalues(RING,'T1',{1});
    ddy     = -atgetfieldvalues(RING,'T1',{3});
    ddphi   = -asin(atgetfieldvalues(RING,'R1',{1,3})); % angle sign to be verified MA 20240130


    % Get indices for the magnet types
    % NB! May be more efficient methods to get these but the cost only
    % applies if plotting
    di = cellfun(@(x) strcmpi(atguessclass(x),'Bend') && x.BendingAngle > 0,RING);
    ri = cellfun(@(x) strcmpi(atguessclass(x),'Bend') && x.BendingAngle < 0,RING);
    qi = cellfun(@(x) strcmpi(atguessclass(x),'Quadrupole'),RING);
    si = cellfun(@(x) strcmpi(atguessclass(x),'Sextupole'),RING);
    oi = cellfun(@(x) strcmpi(atguessclass(x),'Multipole'),RING);
    bpmi = atgetcells(RING,'FamName','BPM');

    figure(2); clf
    subplot(3,1,1); hold on
    plot(spos,ddx,'.'); box on; grid on
    plot(spos(di),ddx(di),'o','color',dcol,'MarkerFaceColor',dcol)
    plot(spos(qi),ddx(qi),'o','color',qcol,'MarkerFaceColor',qcol)
    plot(spos(ri),ddx(ri),'o','color',rcol,'MarkerFaceColor',rcol)
    plot(spos(si),ddx(si),'o','color',scol,'MarkerFaceColor',scol)
    plot(spos(oi),ddx(oi),'o','color',ocol,'MarkerFaceColor',ocol)
    plot(spos(bpmi),ddx(bpmi),'sq','color','k','MarkerFaceColor','k')
    

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
    %     for i = 1:length(qi)
    %         lq(i)  = RING{qi(i)}.Length;
    %         kq(i)  = RING{qi(i)}.PolynomB(2);
    %     end
    lq = atgetfieldvalues(RING(qi),'Length');
    kq = atgetfieldvalues(RING(qi),'PolynomB',{1,2});
    sq = sq+lq/2;
    stem(sq, kq,'linewidth',3,'marker','none')

    sr    = spos(ri);
    lr = atgetfieldvalues(RING(ri),'Length');
    kr = atgetfieldvalues(RING(ri),'PolynomB',{1,2});
    sr = sr + lr/2;
    stem(sr, kr,'linewidth',3,'marker','none')
    axis([0 528/5 -7 7]); xlabel('S (m)'); ylabel('K (m-2)')

    sd    = spos(di);
    ld = atgetfieldvalues(RING(di),'Length');
    kd = atgetfieldvalues(RING(di),'PolynomB',{1,2});
    th = atgetfieldvalues(RING(di),'BendingAngle');

    %     for i = 1:length(di)
    %         ld(i)  = RING{di(i)}.Length;
    %         th(i)  = RING{di(i)}.BendingAngle;
    %         kd(i)  = RING{di(i)}.PolynomB(2);
    %     end
    sd = sd + ld/2;
    stem(sd, kd,'linewidth',3,'marker','none')
    Bdip = th./ld * 3/.2998;

    subplot(3,1,2); hold on
    title('sextupole/octupole components')
    ss    = spos(si);
    %     for i = 1:length(si)
    %         ls(i)   = RING{si(i)}.Length;
    %         k2s(i)  = RING{si(i)}.PolynomB(3);
    %     end
    ls = atgetfieldvalues(RING(si),'Length');
    k2s = atgetfieldvalues(RING(si),'PolynomB',{1,3});

    ss = ss + ls/2; yyaxis left;
    stem(ss, k2s,'linewidth',3,'marker','none');
    axis([0 528/5 -300 300]); ylabel('K (m-3)')
    so    = spos(oi);
    %     for i = 1:length(oi)
    %         lo(i)   = RING{oi(i)}.Length;
    %         k3o(i)  = RING{oi(i)}.PolynomB(4);
    %     end
    lo = atgetfieldvalues(RING(oi),'Length');
    k3o = atgetfieldvalues(RING(oi),'PolynomB',{1,4});

    so = so + lo/2; yyaxis right;
    stem(so, k3o,'linewidth',3,'marker','none')
    axis([0 528/5 -3000 3000]); xlabel('S (m)'); ylabel('K (m-4)')

    subplot(3,1,3); hold on
    title('dipole components')

    %     for i=1:numel(sd)
    %         b(i) = bar(sd(i), Bdip(i), 'BarWidth',ld(i),'Facecolor','cyan');
    %     end
    stairs(sd, Bdip,'b');
    axis([0 528/5 -0.1 0.7]); xlabel('S (m)'); ylabel('B (T)')
end


end



% % % % %     neigeny = 110; neigenx = 160;
% % % % %     [OCS, OCS0, V, S, ErrorFlagx]  = setorbitMA({GoalOrbitX}, {getx('physics','struct',BPMdevlist)},{getsp('HCM','physics','struct')},2,neigenx,'CorrectorGain', 0.85);
% % % % %     [OCS, OCS0, V, S, ErrorFlagy]  = setorbitMA({GoalOrbitY}, {gety('physics','struct',BPMdevlist)},{getsp('VCM','physics','struct')},2,neigeny,'CorrectorGain', 0.85);
% % % % % function v = randn_t(a,b, trunc)
% % % % % outlier = ones(a,1);
% % % % % v = randn(a,b);
% % % % % while sum(outlier)>0;
% % % % %     v(find(outlier>0)) = randn(length(find(outlier>0)),b);
% % % % %     outlier = abs(v)>trunc; % trunc-sigma truncation
% % % % % end
% % % % % end

