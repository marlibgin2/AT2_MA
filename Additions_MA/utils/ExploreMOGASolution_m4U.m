function Rout = ExploreMOGASolution_m4U(MOGAsol, SolN, Xout, Fvalout)
if nargin<1
    MOGAsol = 'MOGA_03112023_070000.mat';
    SolN=3;
    fitchro=0;
end
if nargin<3
load (MOGAsol)
else
    X = Xout;
    fval = Fvalout;
end

figure(11); clf; 
x0=10;
y0=10;
width=1200;
height=400;
set(gcf,'position',[x0,y0,width,height],'color','w')

%
% choose a solution
%

%SolN = 242; 
%figure(11); clf; 
subplot(1,3,1)
plot(fval(:,1)*1e12,fval(:,2)*1e6,'o'); hold on; grid on
XL = get(gca, 'XLim'); YL = get(gca, 'YLim');
text(XL(1)+(XL(2)-XL(1))/4,YL(2)-abs(YL(2)-YL(1))/6,['emix = ' num2str(fval(SolN)*1e12,4) ' (pm.rad)'],'color','r')
text(XL(1)+(XL(2)-XL(1))/4,YL(2)-abs(YL(2)-YL(1))/6-abs(YL(2)-YL(1))/10,['Sol. ' num2str(SolN)],'color','r')
%text(fval(SolN,1)*1e12-1,fval(SolN,2)*1e6+5,['emix = ' num2str(fval(SolN)*1e12,4) ' (pm.rad)'],'color','r')
%text(fval(SolN,1)*1e12-0.5,fval(SolN,2)*1e6-5,['Sol. ' num2str(SolN)],'color','r')
plot(fval(SolN,1)*1e12,fval(SolN,2)*1e6,'ro','markerfacecolor','r')
xlabel('emix (pm.rad)'); ylabel('-AreaDA (mm^2)')
title('Pareto Front')
set(gca,'Position',[0.07 0.12 0.25 0.75])
gcaPF = gca;
%
% recalculate and plot the DA
%
%load('medmax_7BA_2_1_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
%%% load('maxivu_B_blank_AT2.mat','RINGb'); rrr=RINGb;

%RI  = m4U_240114_b01_02_03_02__grd_segmented;
RI  = m4U_240314_b01_02_04_03;    
rrr  = alter_m4U_lattice(X(SolN,:), RI, 'B'); clear RI; 
s = findspos(rrr,1:length(rrr)+1); se = s(end);
% guess the periodicity
P = round(528/se); 
rrr = atfitchrom(rrr,[1/P 1/P],'S1','S2'); 

r0     = 3e-3;
nsteps = 20; %61; 
nturns1 = 250;
dp     = 0.0; 
res    = 0.25e-3;
alpha  = 1.1; 
DA     = modelDA_sim_par(rrr, r0, nsteps, nturns1*P, dp, res, alpha); %%with parfor

% r0     = 3e-3;
% nsteps = 31; %61; 
% nturns = 250;
% dp     = 0.0; 
% res    = 0.25e-3;
% alpha  = 1.1; 
% DA     = modelDA_sim_par(rrr, r0, nsteps, nturns*P, dp, res, alpha);
SumDA  = -calcDA_Area(DA)

figure(11); hold on
subplot(1,3,2)
plot(DA(:,1),DA(:,2),'b-','linewidth',3); grid on; hold on


nturns2 = 2000;
DA_2000t     = modelDA_sim_par(rrr, r0, nsteps, nturns2*P, dp, res, alpha);
SumDA_2000t  = -calcDA_Area(DA_2000t);
plot(DA_2000t(:,1),DA_2000t(:,2),'r:','linewidth',3); 
legend([num2str(nturns1) 't'],[num2str(nturns2) 't'])
xlabel('X (m)'); ylabel('Y (m)'); axis([-0.015 0.015 0 10e-3])
XL = get(gca, 'XLim'); YL = get(gca, 'YLim');
text(XL(1)+(XL(2)-XL(1))/5,YL(2)-abs(YL(2)-YL(1))/5,['sumDA = ' num2str(-SumDA*1e6,8) ' (mm^2)'],'color','b');

title('Dynamic Aperture')
set(gca,'Position',[0.39 0.12 0.25 0.75])

text(-0.015,0.011,[MOGAsol(1:end-4) '_S' num2str(SolN)],'interpreter','none')
gcaDA = gca;

% -------------------------
% plot the Twiss Parameters
% -------------------------

figure(11); hold on
subplot(1,3,3); 
CURVE = atplot(rrr); % [penalty Rout]=fun_medmax_match_AT2(X(SolN,:),rrr, 1, 1);
set(gca,'Position',[0.69 0.12 0.25 0.75])
gcaTWI = gca; 

Rout = rrr; 

saveas(gcf,[MOGAsol(1:end-4) '_Sol' num2str(SolN) '.png'])
saveas(gcf,[MOGAsol(1:end-4) '_Sol' num2str(SolN) '.fig'])
end

function DA_Area = calcDA_Area(DA)
a = []; 
A = [0, 0];
for i = 1:length(DA)-1
    B = DA(i,  :);
    C = DA(i+1,:);
    a(i) = (A(1) * (B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))/2;
end
DA_Area = sum(a);% * W;% * AspRa;
end
