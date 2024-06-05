function [maxAmplitude, IE, Npart0, Npart] = ...
    at_m4U_Injection(Rin, Imik, theta_dk, x0, x0p, N_of_kicks, Nturn, beam_bunch_number,igen)% deltasqfo)

graphic = 1;
MIK_ON = 1;
DK_ON  = 1;

% -------------------
% upload the AT model
% -------------------

maxAmplitude =  0;
IE = 0;

DK                = set_DK_parameters();
MIK               = set_MIK_parameters();

[SB, IB, Rin] = set_machine_parameters(Rin); % input are the ref bunch n. for SB and IB
S = SB.S;

% ----------------------------------------------------------
% wrap it up in 6D vectors (sigmas of the beams in ph-space)
% ----------------------------------------------------------
% S_inj      = [sx_linac sxp_linac 0 0 sE_linac 0]';
fINJ = 1; fSTO = 1;
S_inj      = fINJ * [IB.sx IB.sxp IB.sy IB.syp IB.sE 0]';
S_sto      = fSTO * [SB.sx SB.sxp 0     0      SB.sE 0]';

% -----------------------------------
% load MIK intensity curve (kick map)
% -----------------------------------
if MIK_ON == 1
    %%%MIK_map = load('MIK_KM.mat');
    MIK_map = load('MIK_KM_EXTENDED_XY.mat'); % analitically extended for larger excursions
    % ------------------------------------------------
    % new, 30052024: ensure that the MIK map is loaded
    % ------------------------------------------------
    % brho =  -3.3356*1.5;
    nx = size(MIK_map.hKM); nx=nx(2);%21
    ny = size(MIK_map.vKM); ny=ny(1);%7
    % ----------------------------------------------
    % to redefine the EPU kick-map
    % ----------------------------------------------
    limx = 0.001*floor(nx/2); limy = 0.002*floor(ny/2);
    [XMIK, YMIK] = meshgrid(-limx:0.001:limx,limy:-0.002:-limy);
    
    PXMIK = -1e-3*MIK_map.hKM; % kick strength to be reviewed
    PYMIK =     0*MIK_map.vKM; % "    "        "  "  " 
    MIKindex = findcells(Rin,'FamName','MIK'); % identify MIK index
     
      Rin{MIKindex}.NumX  = nx;
      Rin{MIKindex}.NumY  = ny;
     
      Rin{MIKindex}.XGrid = XMIK;
      Rin{MIKindex}.YGrid = YMIK;
   
    Rin{MIKindex}.PxGrid = PXMIK;
    Rin{MIKindex}.PyGrid = PYMIK;
end

%
% introduce a time series, with a half-sine pulsed function peaking at injection (T1)
% and then silent for T --> oo
%

% ------------------
% general parametersset_vac
% ------------------
%%% Nturn = 7; %7;
Npart = 30; Npart0=Npart;
Npart_sto = 2;

if nargin <=       3.55 % peak current in kA
    theta_dk     = 2.223365;
    x0       = -13.5e-3;
    x0p      = 0;
end

AmpliFac = Imik / 3.55;
% AmpliFac = 0;
R0_inj   = [ x0  x0p   0 0   0 0]';
R0_sto   = [ 0   0     0 0   0 0]';
silence  = 1; % 0 = silence the random component
X0_inj   = silence*randn(6,Npart).*repmat(S_inj,1,Npart)+repmat(R0_inj,1,Npart); % injected beam
load('InjTest_30part.mat')
%load(['BKUP_' num2str(igen) '.mat'])
%X0_sto   = silence*randn(6,Npart_sto).*repmat(S_sto,1,Npart_sto)+repmat(R0_sto,1,Npart_sto);
%%%save(['BKUP_' num2str(igen) '.mat'],'X0_inj')
%%% save(['InjTest_30part.mat'],'X0_inj')
% -----------------------------------------
% define the Vacuum Chamber / BPM positions
% -----------------------------------------
%[xlim_m, xlim_p,  ylim_m, ylim_p] = set_vacuum_chamber(1000, SB);
% sbpm = findspos(Rin,bpmindex);

for AF = AmpliFac

    if graphic == 1
        figure(11); clf; hold on; subplot(2,4,[1 3]); axis([0 528 -0.03 0.03]); set(gcf,'Position',[0 0 1900 700])
        plot([1.252 1.252],[-0.03 0.03],'k--'); text (1.3,  0.025, 'IP'); hold on
        MIKpos = findspos(Rin,MIKindex);
        %plot([28.3 28.3],[-0.03 0.03],'m--'); text (28.8, 0.025, 'MIK','color','m'); hold on
        plot([MIKpos MIKpos],[-0.03 0.03],'m--'); text (MIKpos+0.3, 0.025, 'MIK','color','m'); hold on
        
        plot([5.04 5.04],[-0.03 0.03],'r--'); text (5.55, 0.025, 'DK','color','r'); hold on
        plot([0 528],[-0.004 -0.004],'k--'); hold on
        plot([0 528],[ 0.004  0.004],'k--'); hold on
        xlabel('S (m)'); ylabel('X,Y (m)')
        drawseptum()
        %     stairs(S, xlim_p,'color',[0.5 0.5 0.5],'linewidth',4)
        %     stairs(S, xlim_m,'color',[0.5 0.5 0.5],'linewidth',4)
        axis([0 528 -30e-3 30e-3])
    end

    nt = 0;
    MIKindex = findcells(Rin,'FamName','MIK'); % identify MIK index
    DKindex  = findcells(Rin,'FamName','PMH3'); % identify DK  index

    maxAmplitude     = -1;
    injection_length = 2:300;  % elements corresponding to the initial 2.5m
    for iel = injection_length % corresponds to 2.5m from mid-inj straight
        if isfield(Rin{iel}, 'EApertures' )
            Rin{iel}.EApertures = [0.022, 0.022];
        end
    end

    for t = 0: Nturn
        nt=nt+1;
        if nt==2
            % -----------------------------------------------------
            % restore the small aperture at injection (to mimic the
            % injection out of the septum)
            % -----------------------------------------------------
            for iel = injection_length % corresponds to 2.5m from mid-inj straight
                if isfield(Rin{iel}, 'EApertures' )
                    Rin{iel}.EApertures = [0.011, 0.011];
                    if iel>=4 && iel<=6
                        Rin{iel}.EApertures = [0.01, 0.01];
                    end
                end
            end

        end
        dt = 528 / 299792458 / 176; % time separation between buckets

        if DK_ON == 1
            Rin{DKindex}.KickAngle = DK.func(beam_bunch_number, t, dt) * theta_dk * [1 0];
        end
        if MIK_ON == 1
            Rin{MIKindex}.PxGrid =  MIK.func(beam_bunch_number, t, dt) * AmpliFac * Rin{MIKindex}.PxGrid;
            Rin{MIKindex}.PyGrid =  MIK.func(beam_bunch_number, t, dt) * AmpliFac * Rin{MIKindex}.PyGrid;
        end


        [Rfin, loss, ~]     = linepass(Rin, X0_inj, 1:length(Rin)+1);
        try 
            A = reshape(Rfin,6,Npart,length(Rin)+1); % extract the Npart coordinates
        catch
            disp('bad solution ... set IE to 0')
            IE = 0; Npart = 0; maxAmplitude = 0; 
            return
        end
            B = A(1,:,:); X  = squeeze(B); % permute(B,[3,1,2]);
        B = A(2,:,:); Xp = squeeze(B); % permute(B,[3,1,2]);
        B = A(3,:,:); Y  = squeeze(B); % permute(B,[3,1,2]);
        B = A(4,:,:); Yp = squeeze(B); % permute(B,[3,1,2]);
        B = A(5,:,:); dP = squeeze(B); % permute(B,[3,1,2]);
        B = A(6,:,:); cT = squeeze(B); % permute(B,[3,1,2]);

        if graphic == 1
            figure(11); hold on;
            if t<1
                colore = {'c', 'm', 'co', 'mo'};
            else
                colore = {'b', 'r', 'bo', 'ro'};
            end
            subplot(2,4,[1 3]); plot(S,X,'color',colore{1},'linewidth',1)
            subplot(2,4,[1 3]); plot(S,Y,'color',colore{2},'linewidth',1)

            subplot(2,4,4); plot(X(:,7),Xp(:,7),colore{3},'markersize',1,'markerfacecolor',colore{1}); axis([-0.02 0.02 -3.e-3 3.e-3]); hold on
            subplot(2,4,4); plot(Y(:,7),Yp(:,7),colore{4},'markersize',1,'markerfacecolor',colore{2}); axis([-0.02 0.02 -3.e-3 3.e-3]); hold on
            xlabel('X,Y (m)'); ylabel("X^{'}, Y^{'}(rad)" )
            title('T-phase space @IP')

            subplot(2,4,8); plot(cT(:,7),dP(:,7),colore{3},'markersize',1,'markerfacecolor',colore{1}); axis([-5e-2 5e-2 -10.e-3 10.e-3]); hold on
            ylabel('\deltaP (%)'); xlabel("cT (m)" )
            title('L-phase space @IP')

        end
        % ---------------------------------------------------------------
        % find and eliminate tracks lost in the VC (defined as EAperture)
        % ---------------------------------------------------------------
        X(loss==1,:)=[];
        Xp(loss==1,:)=[];
        Y(loss==1,:)  = [];
        Yp(loss==1,:) = [];
        dP(loss==1,:) = [];
        cT(loss==1,:) = [];
        Npart   = Npart - sum(loss);
        if Npart <=0 
            Npart = 0;
            IE    = 0;
            maxAmplitude=0;
            return
        end
        IE      = 100*(Npart/Npart0);

        X0_inj = [X(:,end)'; Xp(:,end)'; Y(:,end)'; Yp(:,end)'; dP(:,end)'; cT(:,end)']; % Rfin(:,end);
        if nt > 3
            if max(abs(Rfin(1,:))) > maxAmplitude
                maxAmplitude = max(abs(Rfin(1,:)));
            end
        end

    end
end
if graphic==1
    figure(11);  subplot(2,4,[1 3]);hold on;
    text (-80,0.03, ['N_{turns} = ' num2str(Nturn)])
    text (-80,0.025,['IE = ' num2str(IE,3) ' %'])
    axis([0 528 -0.03 0.03])
end
end


function drawseptum()
rectangle('position',[0,-0.0125, 1.252, 0.0025],'FaceColor',[0.5 0.5 0.5])
end

function [SB, IB, Rin] = set_machine_parameters(Rin)
% ------------------------------------------------------
% cavities and radiation should be controlled in the Rin
% ------------------------------------------------------

ATS = atsummary_ring(Rin); % why atsummary(Rin) fails when the ring is "equipped" with DK/MIK?
SB.c = 299792458;
SB.T = ATS.circumference/SB.c;

[TD, ~, ~] = twissring(Rin,0,1:(length(Rin)+1),'chrom',0.0001);
BETA    = cat(1,TD.beta);
ALFA    = cat(1,TD.alpha);
S       = findspos(Rin,1:length(Rin)+1);

% --------------------------------------------
% beam parameters for the STORED BEAM at start
% --------------------------------------------
SB.S     = S;
SB.emix  = ATS.naturalEmittance;
SB.betax = BETA(1,1);
SB.alfax = ALFA(1,1);
SB.gamax = (1+SB.alfax^2)/SB.betax;

SB.sx     = sqrt(SB.emix * SB.betax);
SB.sxp    = sqrt(SB.emix * SB.gamax);
SB.sE     = ATS.naturalEnergySpread; % 7.4505e-04;

% -------------------------------------
% beam parameters just out of the LINAC
% for the INJECTED BEAM (reference?)
% -------------------------------------
InjBeModel = 'Sara';
switch InjBeModel
    case 'paper'
        IB.emix  = 3.4e-9;
        IB.betax = 21.133;
        IB.alfax = -0.002;
        IB.gamax = (1+IB.alfax^2)/IB.betax;

        IB.sx   = sqrt(IB.emix * IB.betax);
        IB.sxp  = sqrt(IB.emix * IB.gamax);

        % y components (guess 23-9-2021)
        IB.emiy  = 10e-12;
        IB.betay = 2;
        IB.alfay = 0;
        IB.gamay = (1+IB.alfay^2)/IB.betay;

        IB.sy   = sqrt(IB.emiy * IB.betay);
        IB.syp  = sqrt(IB.emiy * IB.gamay);

        IB.sE   = 1e-3;
    case 'Sara'  %240307MultibunchUniformE
        % -----------------------------------------------------------------
        % combination of data from above document + some educated guess ...
        % -----------------------------------------------------------------
        IB.emix  = 5.8449e-9;
        IB.betax = 15.3461;
        IB.alfax = -1.0917;
        IB.gamax = (1+IB.alfax^2)/IB.betax;

        IB.sx   = sqrt(IB.emix * IB.betax);
        IB.sxp  = sqrt(IB.emix * IB.gamax);

        % y components (guess 30-5-2024)
        IB.emiy  = 0.01 * IB.emix;
        IB.betay = 4.8059;
        IB.alfay = -1.0191;
        IB.gamay = (1+IB.alfay^2)/IB.betay;

        IB.sy   = sqrt(IB.emiy * IB.betay);
        IB.syp  = sqrt(IB.emiy * IB.gamay);

        IB.sE   = 3e-3;

end
end

function [xlim_m, xlim_p, ylim_m, ylim_p] = set_vacuum_chamber(Rin, vcsizefactor, SB) %#ok<DEFNU> 
% -----------------------------------------
% define the Vacuum Chamber
% -----------------------------------------
S = SB.S; clear SB

if nargin < 1
    vcsizefactor = 1; %1
end
VC           = dlmread('vacuum_chambers_R3_basic.txt','',1,0); % load the vacuum chamber
VC(:,2:end) = VC(:,2:end) * vcsizefactor;
xlim_m = zeros(size(Rin));
xlim_p = zeros(size(Rin));
ylim_m = zeros(size(Rin));
ylim_p = zeros(size(Rin));
for i = 1:length(Rin) % define the VC limits
    %display(['i= ' num2str(i)])
    [a, ~]= find(VC(:,1)>S(i));
    bb   = a(1);
%     [a, ~]= find(VC(:,1)<=S(i));
%     aa   = a(end);
    Rin{i}.Limits = [VC(bb,2) VC(bb,3) VC(bb,4) VC(bb,5) ]; % [-18e-3 18e-3  -5.e-3 5.e-3];
    xlim_m(i)     = Rin{i}.Limits(1);
    xlim_p(i)     = Rin{i}.Limits(2);
    ylim_m(i)     = Rin{i}.Limits(3);
    ylim_p(i)     = Rin{i}.Limits(4);
end
xlim_m(end+1)=xlim_m(end);
xlim_p(end+1)=xlim_p(end);
ylim_m(end+1)=ylim_m(end);
ylim_p(end+1)=ylim_p(end);

end

function DK = set_DK_parameters
DK.t     = 3e-6; % half sine save time legnth
DK.h     = 176; % n. of buvckets
DK.func  = @kick;
DK.dt    = 528 / 299792458 / DK.h; %bunch interdispance
end

function MIK = set_MIK_parameters
MIK.t     = 3e-6; % half sine save time legnth
MIK.h     = 176; % n. of buvckets
MIK.func  = @kick;
MIK.dt    = 528 / 299792458 / MIK.h; %bunch interdistance
end

function Y = kick(nb,nT,dt) %kick(nb, dt)
% h = 176;
T = 528 / 299792458;
for k = 1:length(nb)
    if nb(k)>=0
        %        y = cos((-nT + nb .* dt./T) .* pi/2 *  2*T /3e-6 -nT);
        y = cos((-nT + nb .* dt./T) .* pi/2);
        y(y<0)=0;
    else
        %       y = cos((nT + nb .* dt./T) .* pi/2 *  2*T /3e-6);
        y = cos((nT + nb .* dt./T) .* pi/2);
        y(y<0)=0;
    end
    y(nT>=3)=0;
    Y = y(k);
end
%    h = 176;
%    T = 528 / 299792458;
%    dt = T / h;
%    y = cos((nT + nb .* dt./T) .* pi/2 *  2*T /3e-6);
%    y(y<0)   = 0;
%    y(nT>=3) = 0;
end