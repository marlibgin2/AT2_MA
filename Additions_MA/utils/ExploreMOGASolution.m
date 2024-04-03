function Rout = ExploreMOGASolution(MOGAsol, SolN, Xout, Fvalout)
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
text(fval(SolN,1)*1e12-1,fval(SolN,2)*1e6+5,['emix = ' num2str(fval(SolN)*1e12,4) ' (pm.rad)'],'color','r')
text(fval(SolN,1)*1e12-0.5,fval(SolN,2)*1e6-5,['Sol. ' num2str(SolN)],'color','r')
plot(fval(SolN,1)*1e12,fval(SolN,2)*1e6,'ro','markerfacecolor','r')
xlabel('emix (pm.rad)'); ylabel('-AreaDA (mm^2)')
title('Pareto Front')
set(gca,'Position',[0.07 0.12 0.25 0.75])
gcaPF = gca;
%
% recalculate and plot the DA
%
%load('medmax_7BA_2_1_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
load('maxivu_B_blank_AT2.mat','RINGb'); rrr=RINGb;

oxxoi = findcells(rrr,'FamName','oxxo');
oxyoi = findcells(rrr,'FamName','oxyo');
oyyoi = findcells(rrr,'FamName','oyyo');
VARi = {oxxoi; oxyoi; oyyoi};

% zero all the octupoles )in force from 24-11-2023)
noOCT=1;
if noOCT==1
for j=1:3
    for i=1:length(VARi{j})
        rrr{VARi{j}(i)}.PolynomB(4) = 0;
    end
end
end
[SumDA DA]= obj_modelDA(X(SolN,:), rrr, 100, 1);
figure(11); hold on
subplot(1,3,2)
plot(DA(:,1),DA(:,2),'b-','linewidth',3); grid on; hold on

[SumDA_2000t DA_2000t]= obj_modelDA(X(SolN,:), rrr, 2000, 1);
plot(DA_2000t(:,1),DA_2000t(:,2),'r:','linewidth',3); 
legend('100t','2000t')
xlabel('X (m)'); ylabel('Y (m)'); axis([-0.015 0.015 0 10e-3])
title('Dynamic Aperture')
set(gca,'Position',[0.39 0.12 0.25 0.75])
text(-0.015,0.011,[MOGAsol(1:end-4)],'interpreter','none')
gcaDA = gca;
%
% plot the Twiss Parameters
%
figure(11); hold on
subplot(1,3,3); 
[penalty Rout]=fun_medmax_match_AT2(X(SolN,:),rrr, 1, 1);
set(gca,'Position',[0.69 0.12 0.25 0.75])
gcaTWI = gca; 

saveas(gcf,[MOGAsol(1:end-4) '_Sol' num2str(SolN) '.png'])
saveas(gcf,[MOGAsol(1:end-4) '_Sol' num2str(SolN) '.fig'])
