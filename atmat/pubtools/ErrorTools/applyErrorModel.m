function varargout = applyErrorModel(varargin)%, OrbitCorrectionFlag)
% SPOIL_THE_LATTICE deployes misalignments for sliced magnets and
% [RINGWE, RING0, MAGe, GIRe, Xbpm, Ybpm] = applyErrorModel(RING,ErrorModel,...)
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
% 2. ErrorModel {cell array of structs, optional} See the error model
%               functions for details.
% 3. Flags:
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
%       'NoMisalignments'
%       'NoFieldErrors'
%       'NoGirderErrors'
%
% OUTPUT
% 1. RINGWE -   {cell array of structs} AT2 lattice output, with new errors
%               applied.
% 2. RING0  -   {cell array of structs} AT2 lattice output, without errors.
%
% See also markSlicedMagnets, getmagnetslices, getMagGroupsFromGirderIndex,
% atguessclass


%% Default settings
DisplayFlag = false; % plot final results
ExpandFlag = true;
RelativeFlag = true;
MisalignmentFlag = true;
FieldErrorFlag = true;
GirderFlag = true;

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
            case {'nomisalignments','nomisalignment'}
                MisalignmentFlag = false;
            case {'nofielderror','nofielderrors'}
                FieldErrorFlag = false;
            case {'nogirdererrors','nogirdererror'}
                GirderFlag = false;
        end
    elseif isstruct(varargin{n})
        % Check whether it's an error model...
        if any(isfield(varargin{n},{'Girder','Magnet'}))
            ERRORMODEL = varargin{n};
            varargin(n) = [];
        end
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






%% Deploy girder errors
if GirderFlag
    for n = 1:numel(index_girderelements)
        % Get girder type, if any. Otherwise just apply the baseline
        if isfield(RING{index_girderelements{n}(1)},'GirderType')
            girderType = {RING{index_girderelements{n}(1)}.GirderType};
        else
            girderType = {'Baseline'};
        end
        iFound = filterByGirder(girderType);



        for p = 1:numel(iFound)
            % Depending on the girder type there may be different error magnitudes
            move_grd(index_girderelements{n}, ERRORMODEL.Girder{iFound(p)}.Random{1}, ERRORMODEL.Girder{iFound(p)}.Systematic{1});
        end
    end
end

    function iFound = filterByGirder(girderType)
        iFound = cellfun(@(x) any(strcmpi(x.ID, girderType)), ERRORMODEL.Girder);
        iFound = find(iFound);
    end


%% Deploy magnet errors

% First, filter away all errors that do not have an ID and warn the user
k = findcells(ERRORMODEL.Magnet,'ID');
if numel(k) < numel(ERRORMODEL.Magnet)
    warning('applyErrorModel:Magnet error model contained error specifications not associated to any magnet type (missing ID). These have been ignored.');
end
ERRORMODEL.Magnet = ERRORMODEL.Magnet(k);
EM_classIDs = getcellstruct(ERRORMODEL.Magnet,'ID',k);

for n = 1:numel(index_magnetslices)
    % Get class
    % NB! May be better to base the guess on integrated strengths, which
    % means a new function.
    atclass = atguessclass(RING{index_magnetslices{n}(1)});

    % Check which errors are applicable
    k = find(strcmpi(atclass,EM_classIDs));

    % Deploy all applicable errors via nested functions to speed up the execution (no need
    % to pass the RING variable to the functions)
    % NB! Possible speed improvement by removing a loop (strcmpi does on
    % internally)
    for i = 1:numel(ERRORMODEL.Magnet(k))

        % For each applicable error
        args_misalignment = cell(1,2);
        args_field = cell(1,2);

        if isfield(ERRORMODEL.Magnet{k(i)},'Systematic')
            for j = 1:numel(ERRORMODEL.Magnet{k(i)}.Systematic)
                if any(isfield(ERRORMODEL.Magnet{k(i)}.Systematic{j},{'Heave','Sway','Surge','Pitch','Yaw','Roll'}))
                    args_misalignment{1} = ERRORMODEL.Magnet{k(i)}.Systematic{j};
                end
                if any(isfield(ERRORMODEL.Magnet{k(i)}.Systematic{j},{'PolynomA','PolynomB'}))
                    args_field{1} = ERRORMODEL.Magnet{k(i)}.Systematic{j};
                end
            end
        end
        if isfield(ERRORMODEL.Magnet{k(i)},'Random')
            for j = 1:numel(ERRORMODEL.Magnet{k(i)}.Random)
                if any(isfield(ERRORMODEL.Magnet{k(i)}.Random{j},{'Heave','Sway','Surge','Pitch','Yaw','Roll'}))
                    args_misalignment{2} = ERRORMODEL.Magnet{k(i)}.Random{j};
                end
                if any(isfield(ERRORMODEL.Magnet{k(i)}.Random{j},{'PolynomA','PolynomB'}))
                    args_field{2} = ERRORMODEL.Magnet{k(i)}.Random{j};
                end
            end
        end

        if MisalignmentFlag && ~all(isempty(args_misalignment))
            move_mag(index_magnetslices{n},args_misalignment{:});
        end

        if FieldErrorFlag && ~all(isempty(args_field))
            applyfielderror(index_magnetslices{n},args_field{:},atclass);
        end

        %         % Check whether this error is a misalignment
        %         if (isfield(ERRORMODEL.Magnet{k(i)},'Systematic') && any(isfield(ERRORMODEL.Magnet{k(i)}.Systematic,{'Heave','Sway','Surge','Pitch','Yaw','Roll'}))) ...
        %                 || (isfield(ERRORMODEL.Magnet{k(i)},'Random') && any(isfield(ERRORMODEL.Magnet{k(i)}.Random,{'Heave','Sway','Surge','Pitch','Yaw','Roll'})))
        %             move_mag(index_magnetslices{n},ERRORMODEL.Magnet{k(i)});
        %         end
        %
        %
        %         % Check whether this error is in a field
        %         if (isfield(ERRORMODEL.Magnet{k(i)},'Systematic') && any(isfield(ERRORMODEL.Magnet{k(i)}.Systematic,{'PolynomA','PolynomB'}))) ...
        %                 || (isfield(ERRORMODEL.Magnet{k(i)},'Random') && any(isfield(ERRORMODEL.Magnet{k(i)}.Random,{'PolynomA','PolynomB'})))
        %             chgrad_mag(index_magnetslices{n},ERRORMODEL.Magnet{k(i)},atclass);
        %         end
    end
end










% ----------------------
% spoiling functions ...
% ----------------------
% <MSj> In the interest of speed I reworked the below into nested
% functions, to avoid the overhead of copying the lattice for each girder
% or magnet that should be moved.
% Further work is likely needed here, as atsetshift and atsettilt function
% calls also leads to a significant amount of copying.

    function move_grd(Gi, eGR, eGS)
        % Input checks
        if ~isfield(eGR,'Sway'), eGR.Sway = 0; end
        if ~isfield(eGS,'Sway'), eGS.Sway = 0; end
        if ~isfield(eGR,'Yaw'), eGR.Yaw = 0; end
        if ~isfield(eGS,'Yaw'), eGS.Yaw = 0; end
        if ~isfield(eGR,'Pitch'), eGR.Pitch = 0; end
        if ~isfield(eGS,'Pitch'), eGS.Pitch = 0; end
        if ~isfield(eGR,'Heave'), eGR.Heave = 0; end
        if ~isfield(eGS,'Heave'), eGS.Heave = 0; end
        if ~isfield(eGR,'Roll'), eGR.Roll = 0; end
        if ~isfield(eGS,'Roll'), eGS.Roll = 0; end

        % -------------------------------------
        % define coordinates of girder elements
        % and pivotal point for yaw/pitch moves
        % -------------------------------------
        s1     = spos(Gi(1));
        s2     = spos(Gi(end));
        sm     = (s1+s2)/2; % pivotal centre

        % --------------------------
        % horizontal plane: sway/yaw
        % --------------------------
        sway   = eGS.Sway + eGR.Sway * randn();
        yaw    = eGS.Yaw + eGR.Yaw * randn();
        dx     = sway + yaw*(spos(Gi)-sm);

        % ---------------------------
        % vertical plane: heave/pitch
        % ---------------------------
        heave  = eGS.Heave + eGR.Heave * randn();
        pitch  = eGS.Pitch + eGR.Pitch * randn();
        dy     = heave + pitch*(spos(Gi)-sm);

        % -------------------------------
        % move in (x,y), apply roll (phi)
        % -------------------------------
        dphi   = eGS.Roll + eGR.Roll * randn();

        % For speed, call atshiftelem and attiltelem functions directly
        % rather than calling atsetshift and atsettilt. This avoids the
        % quite large overhead of copying the entire lattice a large number
        % of times.
        for ii = 1:length(Gi)
            RING{Gi(ii)}=atshiftelem(RING{Gi(ii)},dx(ii),dy(ii),'RelativeShift');
            RING{Gi(ii)}=attiltelem(RING{Gi(ii)},dphi,'RelativeTilt');
        end


    end



    function applyfielderror(mi, Es, Er, class)

        % The input checks are likely better placed in the error model
        %         if ~isfield(EM,'Systematic'), EM.Systematic = struct('PolynomA',[],'PolynomB',[]); end
        %         if ~isfield(EM.Systematic,'PolynomA'), EM.Systematic.PolynomA = zeros(1,4); end
        %         if ~isfield(EM.Systematic,'PolynomB'), EM.Systematic.PolynomB = zeros(1,4); end
        %
        %         if ~isfield(EM,'Random'), EM.Random = struct('PolynomA',[],'PolynomB',[]); end
        %         if ~isfield(EM.Random,'PolynomA'), EM.Random.PolynomA = zeros(1,4); end
        %         if ~isfield(EM.Random,'PolynomB'), EM.Random.PolynomB = zeros(1,4); end
        %

        if isempty(Es.PolynomA) && isempty(Es.PolynomB) && isempty(Er.PolynomA) && isempty(Er.PolynomB) && Er.Scaling == 1 && Es.Scaling == 1, return; end
        if isempty(Es), Es = zeros(1,4); end
        if isempty(Er), Er = zeros(1,4); end

        % Identify the main field component for higher multipole scaling
        switch lower(class)
            case 'quadrupole'
                MainComponentField = 'PolynomB';
                MainComponentIndex = 2;
            case 'sextupole'
                MainComponentField = 'PolynomB';
                MainComponentIndex = 3;
            case 'multipole'
                MainComponentField = 'PolynomB';
                MainComponentIndex = 4;
            case 'corrector'
                MainComponentField = 'PolynomB';
                MainComponentIndex = 1;
            case 'bend'
                MainComponentField = 'BendingAngle';
                MainComponentIndex = 1;
        end

        maxOrder = max([ numel(Er.PolynomB), numel(Er.PolynomA), numel(Es.PolynomB), numel(Es.PolynomA)]);

        % Pad to the max order of error specified
        Er.PolynomB = padZeros(Er.PolynomB,maxOrder);
        Er.PolynomA = padZeros(Er.PolynomA,maxOrder);
        Es.PolynomB = padZeros(Es.PolynomB,maxOrder);
        Es.PolynomA = padZeros(Es.PolynomA,maxOrder);

        function y = padZeros(x,N)
            pad = zeros(1,N); pad(1:numel(x)) = x; y = pad;
        end

        % Generate the error for this magnet
        Et.PolynomB = Es.PolynomB + Er.PolynomB .* trunc_randn(maxOrder,2)';
        Et.PolynomA = Es.PolynomA + Er.PolynomA .* trunc_randn(maxOrder,2)';

        % Apply error
        for ii = 1:numel(mi)
            RING{mi(ii)}.PolynomB = padZeros(RING{mi(ii)}.PolynomB,maxOrder) + RING{mi(ii)}.(MainComponentField)(MainComponentIndex) * Et.PolynomB;
            RING{mi(ii)}.PolynomA = padZeros(RING{mi(ii)}.PolynomA,maxOrder) + RING{mi(ii)}.(MainComponentField)(MainComponentIndex) * Et.PolynomA;
            RING{mi(ii)}.MaxOrder = maxOrder - 1;
        end

    end


% -----------------------
% move individual magnets
% -----------------------

    function move_mag(mi, Es, Er)
        % ------------------------------------------
        % input: mi, magnet index / me: magnet error
        % ------------------------------------------
        % Generate shifts and rolls for this single magnet
        dx     = Es.Sway + Er.Sway * trunc_randn(1,2);
        dy     = Es.Heave + Er.Heave * trunc_randn(1,2);
        dphi   = Es.Roll + Er.Roll * trunc_randn(1,2);

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

