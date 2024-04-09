% macro match dba test lattice beta functions and dispersion using
% quadrupoles.
%
% this macro shows the available functionalities of atmatch. 
% 
% various variable and constraint input constructions are shown

clear all

%load('medmax_7BA_1_1_1_AT2.mat','RING','RING1');RING=RING1; 
%%% load('medmax_7BA_1_1_1_AT2_M00.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M1.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M2.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M3.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M4.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M5.mat','RING','RING1');RING=RING1; 
% load('medmax_7BA_1_1_1_AT2_M6.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M7.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M8.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M9.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M12.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M13.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M14.mat','RING','RING1');RING=RING1; 
%load('medmax_7BA_1_1_1_AT2_M14.mat','RING','RING1');RING=RING1; % ex=250pm
%load('medmax_7BA_1_1_1_AT2_M16.mat','RING','RING1');RING=RING1; % ex=250pm

% restart from start 
% load('medmax_7BA_1_1_1_AT2_M18.mat','RING','RING1');RING=RING1; % ex=310pm DIP/DIPm +1.0%
%load('medmax_7BA_1_1_1_AT2_M19.mat','RING','RING1');RING=RING1; % ex=290pm DIP/DIPm +2.0%
%load('medmax_7BA_1_1_1_AT2_M20.mat','RING','RING1');RING=RING1; % ex=270pm DIP/DIPm +3.0%
%load('medmax_7BA_1_1_1_AT2_M21.mat','RING','RING1');RING=RING1; % ex=265pm DIP/DIPm +4.0%
%load('medmax_7BA_1_1_1_AT2_M22.mat','RING','RING1');RING=RING1; % ex=260pm DIP/DIPm +4.5%
%load('medmax_7BA_1_1_1_AT2_M23.mat','RING','RING1');RING=RING1; % ex=257pm DIP/DIPm +5.0%
%load('medmax_7BA_1_1_1_AT2_M24.mat','RING','RING1');RING=RING1; % ex=257pm DIP/DIPm +5.0%
%load('medmax_7BA_1_1_1_AT2_M25.mat','RING','RING1');RING=RING1; % ex=257pm DIP/DIPm +5.1%
%load('medmax_7BA_1_1_1_AT2_M26.mat','RING','RING1');RING=RING1; % ex=257pm DIP/DIPm +5.2%
%load('medmax_7BA_1_1_1_AT2_M27.mat','RING','RING1');RING=RING1; % ex=255pm DIP/DIPm +5.1%
%load('medmax_7BA_1_1_1_AT2_M28.mat','RING','RING1');RING=RING1; % ex=255pm DIP/DIPm +5.1%
%load('medmax_7BA_1_1_1_AT2_M29.mat','RING','RING1');RING=RING1; % ex=255pm DIP/DIPm +5.2%
%load('medmax_7BA_1_1_1_AT2_M30.mat','RING','RING1');RING=RING1; % ex=255pm DIP/DIPm +5.4%
%load('medmax_7BA_1_1_1_AT2_M31.mat','RING','RING1');RING=RING1; % ex=255pm DIP/DIPm +5.6%
%load('medmax_7BA_1_1_1_AT2_M31.mat','RING','RING1');RING=RING1; % ex=254pm DIP/DIPm +5.6%
%load('medmax_7BA_1_1_1_AT2_M33.mat','RING','RING1');RING=RING1; % ex=254pm DIP/DIPm +5.8%
%load('medmax_7BA_1_1_1_AT2_M34.mat','RING','RING1');RING=RING1; % ex=254pm DIP/DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M35.mat','RING','RING1');RING=RING1; % ex=253pm DIP/DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M36.mat','RING','RING1');RING=RING1; % ex=253pm DIP +6.6% /DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M37.mat','RING','RING1');RING=RING1; % ex=253pm DIP +7.1% /DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M38.mat','RING','RING1');RING=RING1; % ex=252pm DIP +7.1% /DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M39.mat','RING','RING1');RING=RING1; % ex=252pm DIP +7.6% /DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M40.mat','RING','RING1');RING=RING1; % ex=253pm DIP +7.6% /DIPm +6.1%
%load('medmax_7BA_1_1_1_AT2_M41.mat','RING','RING1');RING=RING1; % ex=253pm DIP +7.9% /DIPm +6.1%


load('medmax_7BA_1_1_1_AT2_simple.mat','RING','RING1');RING=RING1; 

addpath(fullfile(pwd,'..'))

%%  VARIABLES

% kDIPmstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','DIPm'),1,2);
% indDIPm  = findcells(RING,'FamName','DIPm'); 
% for i=1:length(indDIPm)
%     RING{indDIPm(i)}.GFactor=1;
% end
% VDIPmK = struct('Indx',{findcells(RING,'FamName','DIPm'), ...
%                     @(RING,GradFac)VaryDIPGrad(RING,GradFac,kDIPmstart,'DIPm')},...
%                     'Parameter',{{'GFactor'},1},...
%             'LowLim',{[0.995]},...
%             'HighLim',{[1.005]}...
%             );

% kDIPstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','DIP'),1,2);
% indDIP  = findcells(RING,'FamName','DIP'); 
% for i=1:length(indDIP)
%     RING{indDIP(i)}.GFactor=1;
% end
% VDIPK = struct('Indx',{findcells(RING,'FamName','DIP'), ...
%                     @(RING,GradFac)VaryDIPGrad(RING,GradFac,kDIPstart,'DIP')},...
%             'Parameter',{{'GFactor'},1},...
%             'LowLim',{[0.995]},...
%             'HighLim',{[1.005]}...
%             );

% kDIPstart0 = getcellstruct(RING,'PolynomB',findcells(RING,'FamName','DIP'),1,2);
% %LS6        = (linspace(6,1200-12+6,1200/12));
% VDIPK2 = struct('Indx',{findcells(RING,'FamName','DIP'), ...
%                     @(RING,kDIP6)VaryDIPGrad2(RING,kDIP6,kDIPstart0,'DIP')},...
%             'Parameter',{{'PolynomB',{1,2}},kDIPstart0(6)},...
%             'LowLim',{[kDIPstart0(6)*1.00]},...
%             'HighLim',{[kDIPstart0(6)*1.00]}...
%             );

indDIP  = findcells(RING,'FamName','DIP'); 
for i=1:length(indDIP)
     RING{indDIP(i)}.GFactor=1;
end
VDIPK3 = struct('Indx',{findcells(RING,'FamName','DIP'), ...
                    @(RING,GF)VaryDIPGrad3(RING,GF,'DIP')},...
            'Parameter',{{'GFactor'},1},...     
            'LowLim',{[0.9]},...
            'HighLim',{[1.1]}...
            );

% kDIPmstart0 = getcellstruct(RING,'PolynomB',findcells(RING,'FamName','DIPm'),1,2);
% %LSm6        = (linspace(6,1200-12+6,1200/12));
% VDIPmK2 = struct('Indx',{findcells(RING,'FamName','DIPm'), ...
%                     @(RING,kDIPm6)VaryDIPGrad2(RING,kDIPm6,kDIPmstart0,'DIPm')},...
%             'Parameter',{{'PolynomB',{1,2}},kDIPmstart0(6)},...
%             'LowLim',{[kDIPmstart0(6)*1.00]},...
%             'HighLim',{[kDIPmstart0(6)*1.00]}...
%             );


kQFMstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QFM'),1,2);
VQFMK = struct('Indx',{findcells(RING,'FamName','QFM'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QFM')},...
            'Parameter',{{'PolynomB',{1,2}},kQFMstart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQFstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QF'),1,2);
VQFK = struct('Indx',{findcells(RING,'FamName','QF'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QF')},...
            'Parameter',{{'PolynomB',{1,2}},kQFMstart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQDENDstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QDEND'),1,2);
VQDENDK = struct('Indx',{findcells(RING,'FamName','QDEND'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QDEND')},...
            'Parameter',{{'PolynomB',{1,2}},kQDENDstart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQFENDstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QFEND'),1,2);
VQFENDK = struct('Indx',{findcells(RING,'FamName','QFEND'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QFEND')},...
            'Parameter',{{'PolynomB',{1,2}},kQFENDstart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );
 
onlymatch=1;
REF=1;
if onlymatch ==1
    %RING   = VaryDIPGrad(RING,1.02,kDIPstart,'DIP'); % increase the DIP gradient
    %RING   = VaryDIPGrad(RING,1.005,kDIPmstart,'DIPm'); % increase the DIP gradient
    % RING   = VaryDIPGrad(RING,1.003,kDIPstart,'DIP'); % increase the DIP gradient
    %RING   = VaryDIPGrad(RING,1.003,kDIPmstart,'DIPm'); % increase the DIP gradient
    VarK   = [VQFENDK, VQDENDK, VQFMK, VQFK, VDIPK3];
    %VarAng = [VABmangle];
else
    VarK   = [VABmK, VTGBmK, VQM0K, VQM1K, VQM1aK, VQM2K, VQM3K, VQM4K, VABK, VTGBK];
    %VarAng = [VABmangle, VABangle ];
end
full=0;
Variab = [VarK]; % fit only the gradients
if full==1
%Variab = [VarK, VarAng] ; % fit gradients and angles
end

%Variab = [Variab, VQM0K, VQM1K, VQM1aK, VQM2K, VQM3K, VQM4K, VABK, VTGBK];

%%  CONSTRAINTS


Ax=struct('Fun',@(RING,~,~)alfx(RING,s1hindx(4)),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[1],...
    'Weight',0.1);

Ay=struct('Fun',@(RING,~,~)alfy(RING,s1hindx(4)),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[],...
    'Weight',1e-3);
 
% -----------
% constraints
% -----------
BSx=struct('Fun',@(RING,~,~)betx_MA(RING,1),...
    'Min',1.5,...  %1.5
    'Max',14.5,...  %1.5
    'RefPoints',[1],...
    'Weight',0.1);
BSy=struct('Fun',@(RING,~,~)bety_MA(RING,1),...
    'Min',1.5,... %5
    'Max',4.5,... %5
    'RefPoints',[1],...
    'Weight',1);
BSxy =struct('Fun',@(RING,~,~)abs(betx_MA(RING,1)-bety_MA(RING,1)),...
    'Min',0.0,... %5
    'Max',0.0,... %5D
    'RefPoints',[1],...
    'Weight',2);
DSx=struct('Fun',@(RING,~,~)dispx2(RING,1),...
    'Min',0.0,...  %1.5
    'Max',0.0,...  %1.5
    'RefPoints',[1],...
    'Weight',0.01);
DPSx=struct('Fun',@(RING,~,~)dispxp2(RING,1),...
    'Min',0.0,...  %1.5
    'Max',0.0,...  %1.5
    'RefPoints',[1],...
    'Weight',0.1);
MBy=struct('Fun',@(RING,~,~)max(bety_MA(RING,50:327)),...
    'Min',0,...
    'Max',28,...
    'RefPoints',[1],...
    'Weight',1);
MDx=struct('Fun',@(RING,~,~)max(dispx2(RING,108:279)),...%50:327)),...
    'Min',0.00,...
    'Max',0.06,...
    'RefPoints',[1],...
    'Weight',1e-3);
mDx=struct('Fun',@(RING,~,~)min(dispx2(RING,52:114)),...
    'Min',0.0025,...
    'Max',0.0045,...
    'RefPoints',[],...
    'Weight',1e-3);
EMIX=struct('Fun',@(RING,~,~)emix(RING),...
    'Min',120.e-12,...
    'Max',330.00e-12,...
    'RefPoints',[1],...
    'Weight',1e-13);
% % % NUX=struct('Fun',@(RING,~,~)nux(RING,20),...
% % %     'Min',64.2,... %64.2
% % %     'Max',64.2,...
% % %     'RefPoints',[],...
% % %     'Weight',0.01*1);
% % % NUY=struct('Fun',@(RING,~,~)nuy(RING,20),...
% % %     'Min',26.34,... %26.34
% % %     'Max',26.34,...
% % %     'RefPoints',[],...
% % %     'Weight',0.01*1);
% % % DNUX_uc = struct('Fun',@(RING,~,~)DnuxUC(RING,20),...
% % %     'Min',0.395,...
% % %     'Max',0.395,...
% % %     'RefPoints',[],...
% % %     'Weight',0.01);
% % % DNUY_uc = struct('Fun',@(RING,~,~)DnuyUC(RING,20),...
% % %     'Min',0.14,...
% % %     'Max',0.14,...
% % %     'RefPoints',[],...
% % %     'Weight',0.01);
% % % 
% % % XIX=struct('Fun',@(RING,~,~)xix(RING,20),...
% % %     'Min',1,...
% % %     'Max',1,...
% % %     'RefPoints',[],...
% % %     'Weight',1);
% % % XIY=struct('Fun',@(RING,~,~)xiy(RING,20),...
% % %     'Min',1,...
% % %     'Max',1,...
% % %     'RefPoints',[],...
% % %     'Weight',1);

% Constr = [ D30 EMIX MB22 NUY NUX XIY XIX A12 A21 B31];%  DD202 DD203 ];%MB22 MD22 DD202 XIX XIY]; %  C41 C42 C51 C52 C61 C62];
% Constr = [ D30 MB22 NUX NUY DD202 DD203 A12 A21 EMIX DNUY_uc ];%MB22 MD22 DD202 XIX XIY]; %  C41 C42 C51 C52 C61 C62];

REF = 1;
if REF == 1
    Constr = [ DSx EMIX ];%MB22 MD22 DD202 XIX XIY]; %  C41 C42 C51 C52 C61 C62];
end


%RING_matched_optconstr=atmatch(RING,Variab,Constr,10^-6,5000,3,@lsqnonlin,'verbose',3,'UseParallel',true);%
RING_matched_optconstr=atmatch(RING,Variab,Constr,10^-5,1000,3,@lsqnonlin,'verbose',3,'UseParallel',false);%
%RING_matched_optconstr=atmatch(RING,Variab,Constr,10^-5,1000,3,'verbose',3);%

figure;
atplot(RING); drawnow
xlim([0 528/20]);
% export_fig('ringdba.pdf','-transparent');
%figure;atplot(RING_matched);% export_fig('ringdba_matched.pdf','-transparent');
figure;
atplot(RING_matched_optconstr); drawnow
xlim([0 528/20])% export_fig('ringdba_matched.pdf','-transparent');