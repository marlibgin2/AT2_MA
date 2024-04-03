% macro match dba test lattice beta functions and dispersion using
% quadrupoles.
%
% this macro shows the available functionalities of atmatch. 
% 
% various variable and constraint input constructions are shown

clear all
load('mTME19_match.mat','RING');

addpath(fullfile(pwd,'..'))

%%  VARIABLES

kABmstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','ABm'),1,2);
V1 = struct('Indx',{findcells(RING,'FamName','ABm'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'ABm')},...
            'Parameter',{{'PolynomB',{1,2}},kABmstart(1)},...
            'LowLim',{[]},...
            'HighLim',{[]}...
            );
kTGBmstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','TGBm'),1,2);
V2 = struct('Indx',{findcells(RING,'FamName','TGBm'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'TGBm')},...
            'Parameter',{{'PolynomB',{1,2}},kTGBmstart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQM1start=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QM1'),1,2);
V3 = struct('Indx',{findcells(RING,'FamName','QM1'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QM1')},...
            'Parameter',{{'PolynomB',{1,2}},kQM1start(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQM1astart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QM1a'),1,2);
V3a = struct('Indx',{findcells(RING,'FamName','QM1a'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QM1a')},...
            'Parameter',{{'PolynomB',{1,2}},kQM1astart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQM2start=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QM2'),1,2);
V4 = struct('Indx',{findcells(RING,'FamName','QM2'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QM2')},...
            'Parameter',{{'PolynomB',{1,2}},kQM2start(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQM3start=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QM3'),1,2);
V5 = struct('Indx',{findcells(RING,'FamName','QM3'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QM3')},...
            'Parameter',{{'PolynomB',{1,2}},kQM3start(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kQM4start=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QM4'),1,2);
V6 = struct('Indx',{findcells(RING,'FamName','QM4'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'QM4')},...
            'Parameter',{{'PolynomB',{1,2}},kQM4start(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );

kABstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','AB'),1,2);
V7 = struct('Indx',{findcells(RING,'FamName','AB'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'AB')},...
            'Parameter',{{'PolynomB',{1,2}},kABstart(1)},...
            'LowLim',{[]},...
            'HighLim',{[]}...
            );
kTGBstart=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','TGB'),1,2);
V8 = struct('Indx',{findcells(RING,'FamName','TGB'), ...
                    @(RING,K1Val)VaryQuadFam(RING,K1Val,'TGB')},...
            'Parameter',{{'PolynomB',{1,2}},kTGBstart(1)},...
            'LowLim',{[],[]},...
            'HighLim',{[],[]}...
            );


%Variab1=atVariableBuilder(RING,{'QM1','QM2'},{{'PolynomB',{1,2}},{'PolynomB',{1,2}}});
%k1start=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QDM'),1,2);

%Variab2=atVariableBuilder(RING,{'QM3'},{{'PolynomB',{1,2}}});
%k1start=getcellstruct(RING,'PolynomB',findcells(RING,'FamName','QDM'),1,2);

%Variab3=atVariableBuilder(RING,{'ABm'},{{'PolynomB',{1,2}}});

%Variab4=atVariableBuilder(RING,{'TGBm'},{{'PolynomB',{1,2}}});

% Variab5=atVariableBuilder(RING,{'AB'},{{'PolynomB',{1,2}}});
% 
% Variab6=atVariableBuilder(RING,{'TGB'},{{'PolynomB',{1,2}}});
% 


% Variab2=struct('Indx',{findcells(RING,'FamName','QFM'),@(RING,K1Val)VaryQuadFam(RING,K1Val,'QDM')},...
%     'LowLim',{[],[]},...
%     'HighLim',{[],[]},...
%     'Parameter',{{'PolynomB',{1,2}},k1start(1)}...
%     );

%Variab=[Variab1,Variab2, Variab3, Variab4, Variab5, Variab6];
%Variab=[V1, V2, V3, V3a, V4, V5, V6, V7, V8];
Variab=[V1, V2, V3, V3a,  V4, V5, V6, V7, V8];


%%  CONSTRAINTS
lgbhindx = findcells(RING,'FamName','LGBh');
s1hindx  = findcells(RING,'FamName','S1h');
qm3indx  = findcells(RING,'FamName','QM3');

% original
% C11=struct('Fun',@(RING,~,~)betx(RING,s1hindx(4)),...
%     'Min',4.55006,...
%     'Max',4.55006,...
%     'RefPoints',[],...
%     'Weight',1);

B110=struct('Fun',@(RING,~,~)betx(RING,s1hindx(4)),...
    'Min',4,...
    'Max',4,...
    'RefPoints',[1],...
    'Weight',1);
B210=struct('Fun',@(RING,~,~)bety(RING,s1hindx(4)),...
    'Min',1.1,...
    'Max',1.1,...
    'RefPoints',[1],...
    'Weight',1);

Ds1h4=struct('Fun',@(RING,~,~)dispx2(RING,s1hindx(4)),...
    'Min',0.045,...%0.033
    'Max',0.045,...%0.043
    'RefPoints',[1],...
    'Weight',1e-3);
DP110=struct('Fun',@(RING,~,~)dispxp2(RING,s1hindx(4)),...
    'Min',0,...
    'Max',0,...
    'RefPoints',[1],...
    'Weight',1);

B111=struct('Fun',@(RING,~,~)betx(RING,s1hindx(6)),...
    'Min',4,...
    'Max',4,...
    'RefPoints',[1],...
    'Weight',1);
Ds1h6=struct('Fun',@(RING,~,~)dispx2(RING,s1hindx(6)),...
    'Min',0.033,...%0.033
    'Max',0.033,...%0.043
    'RefPoints',[1],...
    'Weight',1e-3);
DP111=struct('Fun',@(RING,~,~)dispxp2(RING,s1hindx(6)),...
    'Min',0,...
    'Max',0,...
    'RefPoints',[1],...
    'Weight',1);

B112=struct('Fun',@(RING,~,~)betx(RING,s1hindx(8)),...
    'Min',4,...
    'Max',4,...
    'RefPoints',[1],...
    'Weight',1);
D112=struct('Fun',@(RING,~,~)dispx2(RING,s1hindx(8)),...
    'Min',0.033,...
    'Max',0.043,...
    'RefPoints',[1],...
    'Weight',1);
DP112=struct('Fun',@(RING,~,~)dispxp2(RING,s1hindx(8)),...
    'Min',0,...
    'Max',0,...
    'RefPoints',[1],...
    'Weight',1);
A12=struct('Fun',@(RING,~,~)alfx(RING,s1hindx(4)),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[1],...
    'Weight',0.1);

DB110=struct('Fun',@(RING,~,~)abs(bety(RING,lgbhindx(7)) - bety(RING,lgbhindx(5))),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[1],...
    'Weight',1e-2);

D200=struct('Fun',@(RING,~,~)dispx2(RING,lgbhindx(1)),...
    'Min',0.0110,...
    'Max',0.0114,...
    'RefPoints',[1],...
    'Weight',1);
DD200=struct('Fun',@(RING,~,~)abs(dispx2(RING,lgbhindx(5))-dispx2(RING,lgbhindx(3))),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[1],...
    'Weight',3e-3);

DP200=struct('Fun',@(RING,~,~)dispxp2(RING,lgbhindx([1])),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[],...
    'Weight',5e-3);
D201=struct('Fun',@(RING,~,~)dispx2(RING,lgbhindx([3])),...
    'Min',0.0064,...
    'Max',0.007,...
    'RefPoints',[],...
    'Weight',1);
DP201=struct('Fun',@(RING,~,~)dispxp2(RING,lgbhindx([3])),...
    'Min',0.0674,...
    'Max',0.08,...
    'RefPoints',[],...
    'Weight',1);
DD201=struct('Fun',@(RING,~,~)abs(dispx2(RING,lgbhindx(7))-dispx2(RING,lgbhindx(5))),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[1],...
    'Weight',1e-4);
DD202=struct('Fun',@(RING,~,~)abs(dispx2(RING,s1hindx(5))-dispx2(RING,s1hindx(3))),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[1],...
    'Weight',1e-4);

% original
% C21=struct('Fun',@(RING,~,~)bety(RING,s1hindx(4)),...
%     'Min',1.64032,...
%     'Max',1.64032,...
%     'RefPoints',[],...
%     'Weight',1);
B21=struct('Fun',@(RING,~,~)bety(RING,s1hindx(4)),...
    'Min',1.1,...
    'Max',1.1,...
    'RefPoints',[],...
    'Weight',1);
A21=struct('Fun',@(RING,~,~)alfy(RING,s1hindx(4)),...
    'Min',0.0,...
    'Max',0.0,...
    'RefPoints',[],...
    'Weight',1e-3);
MB22=struct('Fun',@(RING,~,~)max(bety(RING,1:166)),...
    'Min',20,...
    'Max',27,...
    'RefPoints',[],...
    'Weight',1);
MD22=struct('Fun',@(RING,~,~)max(dispx2(RING,1:166)),...
    'Min',0.0356,...
    'Max',0.0356,...
    'RefPoints',[],...
    'Weight',5e-3);
mD22=struct('Fun',@(RING,~,~)min(dispx2(RING,1:166)),...
    'Min',0.0065,...
    'Max',0.0065,...
    'RefPoints',[],...
    'Weight',1e-3);

 D30=struct('Fun',@(RING,~,~)dispx2(RING,1),...
     'Min',0.000,...
     'Max',0.000,...
     'RefPoints',[1],...
     'Weight',1e-4);

 D40=struct('Fun',@(RING,~,~)dispx2(RING,166),...
     'Min',1e-7,...
     'Max',1e-7,...
     'RefPoints',[1],...
     'Weight',1);
 
 % C30=struct('Fun',@(RING,~,~)dispx(RING,1),...
 %     'Min',1e-2,...
 %     'Max',1e-2,...
 %     'RefPoints',[1],...
 %     'Weight',1); 
B31=struct('Fun',@(RING,~,~)betx(RING,1),...
    'Min',3,...  %1.5
    'Max',3,...  %1.5
    'RefPoints',[1],...
    'Weight',1);
B32=struct('Fun',@(RING,~,~)bety(RING,1),...
    'Min',3,... %5
    'Max',3,... %5
    'RefPoints',[1],...
    'Weight',1);
Bxy30 =struct('Fun',@(RING,~,~)abs(betx(RING,1)-bety(RING,1)),...
    'Min',0.0,... %5
    'Max',0.0,... %5
    'RefPoints',[1],...
    'Weight',1e-2);

C41=struct('Fun',@(RING,~,~)betx(RING,lgbhindx(9)),...
    'Min',.48579,...
    'Max',.48579,...
    'RefPoints',[1],...
    'Weight',1);
C42=struct('Fun',@(RING,~,~)bety(RING,lgbhindx(9)),...
    'Min',4.8018,...
    'Max',4.8018,...
    'RefPoints',[1],...
    'Weight',1);
C51=struct('Fun',@(RING,~,~)betx(RING,lgbhindx(7)),...
    'Min',.48579,...
    'Max',.48579,...
    'RefPoints',[],...
    'Weight',1);
C52=struct('Fun',@(RING,~,~)bety(RING,lgbhindx(7)),...
    'Min',4.8018,...
    'Max',4.8018,...
    'RefPoints',[1],...
    'Weight',1);
C61=struct('Fun',@(RING,~,~)betx(RING,lgbhindx(5)),...
    'Min',.48579,...
    'Max',.48579,...
    'RefPoints',[1],...
    'Weight',1);
C62=struct('Fun',@(RING,~,~)bety(RING,lgbhindx(5)),...
    'Min',4.8018,...
    'Max',4.8018,...
    'RefPoints',[1],...
    'Weight',1);
EMIX=struct('Fun',@(RING,~,~)emix(RING),...
    'Min',33.e-12,...
    'Max',33.e-12,...
    'RefPoints',[],...
    'Weight',1e-13*100);
NUX=struct('Fun',@(RING,~,~)nux(RING,20),...
    'Min',64.45,...
    'Max',64.45,...
    'RefPoints',[],...
    'Weight',0.01*100);
NUY=struct('Fun',@(RING,~,~)nuy(RING,20),...
    'Min',24.28,...
    'Max',24.28,...
    'RefPoints',[],...
    'Weight',0.01*100);


%Constr = [C11, C12, C21, C22];
Constr = [ B110 B210 DB110 Ds1h4 DP110 B111 DP111 ...
           B112 D112 DP112 A12 A21 D30 DP200 DD200 ...
           DP201 DD201 DD202 MB22 Bxy30 ...
           EMIX NUX NUY]; %  C41 C42 C51 C52 C61 C62];

%Constr = [C11]
%Constr = [C112 ]
if 1==0

Constr1=struct('Fun',@(RING,~,~)dispx(RING,1),...
    'Min',0e-3,...
    'Max',0e-3,...
    'RefPoints',[1],...
    'Weight',1);
disp('Horizontal dispersion at straigth section= 6e-3')

Constr2=struct('Fun',@(RING,~,~)betx(RING,lgbhindx(3)),...
    'Min',.48579,...
    'Max',.48579,...
    'RefPoints',[lgbhindx(3)],...
    'Weight',1);
disp('Horizontal beta at LGBh = 0.48')
Constr3=struct('Fun',@(RING,~,~)alfx(RING,lgbhindx(3)),...
    'Min',0,...
    'Max',0,...
    'RefPoints',[lgbhindx(3)],...
    'Weight',1);
disp('Horizontal alfa at LGBh = 0.00')

Constr4=struct('Fun',@(RING,~,~)bety(RING,lgbhindx(3)),...
    'Min',4.8018,...
    'Max',4.8018,...
    'RefPoints',[lgbhindx(3)],...
    'Weight',1);
disp('Vertical beta at LGBh = 4.8')
Constr5=struct('Fun',@(RING,~,~)alfy(RING,lgbhindx(3)),...
    'Min',0.00,...
    'Max',0.00,...
    'RefPoints',[lgbhindx(3)],...
    'Weight',1);
disp('Vertical alfa at LGBh = 0.00')

Constr6=struct('Fun',@(RING,~,~)dispx(RING,lgbhindx(3)),...
    'Min',0.00622479,...
    'Max',0.00622479,...
    'RefPoints',[lgbhindx(3)],...
    'Weight',1);
disp('Horizontal dispersion at LGBh = 0.00622479')

%qm1indx=findcells(RING,'FamName','QM1');


% Constr3=struct('Fun',{@(RING,~,~)bety(RING,qm1indx(1)),@(~,ld,~)mux(ld)},...
%     'Min',{4,2.5},...
%     'Max',{4,2.5},...
%     'RefPoints',{[],[1:length(RING)+1]},...
%     'Weight',{1,1});
% disp('Vertical beta at QM1= 4')
% disp('Horizontal phase advance = 2.8')

%Constr=[Constr1,Constr2,Constr3,Constr4,Constr5,Constr6];

%% MATCHING
 disp('wait few iterations')

 RING_matched=atmatch(RING,Variab,Constr,10^-20,1000,3,@lsqnonlin);
end

%return
% c1=atlinconstraint(lgbhindx(3),...
%     {{'beta',{1}},{'beta',{2}}},...
%     [0.48579,4.8018],...
%     [0.48579,4.8018],...
%     [1 1]);
% 
% c2=atlinconstraint(lgbhindx(3),...
%     {{'alpha',{1}},{'alpha',{2}}},...
%     [0.0,0.0],...
%     [0.0,0.0],...
%     [1e3 1e3]);
% 
% c3=atlinconstraint(lgbhindx(3),...
%     {{'Dispersion',{1}}},...
%      0.00622479,...
%      0.00622479,...
%      1 );
% 
% c4=atlinconstraint(1,...
%     {{'Dispersion',{1}}},...
%      0 ,...
%      0 ,...
%      1e6 );
% 
% c=[c1,c2,c3,c4];

RING_matched_optconstr=atmatch(RING,Variab,Constr,10^-5,5000,3,@lsqnonlin);% ,'UseParallel',true);%
% RING_matched_optconstr=atmatch(RING,Variab,Constr,10^-8,1000,3);%

figure;atplot(RING);% export_fig('ringdba.pdf','-transparent');
%figure;atplot(RING_matched);% export_fig('ringdba_matched.pdf','-transparent');
figure;atplot(RING_matched_optconstr);% export_fig('ringdba_matched.pdf','-transparent');