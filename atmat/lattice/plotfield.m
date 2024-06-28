function plotfield(varargin)
% plots the magnetic field with their longitudinal position 
%   
%% Inputs
% Mandatory argument
% FG : structure generated by calcMagneticFields
%
%% Usage examples
% plotfield(FG);

%% History
% SJ  2024/06/27: first version
% PFT 2024/06/28: added grid lines and title

FG       = getargs(varargin,[]);
desc=FG.desc;
Spos = FG.Spos;
field= FG.Field;
GradQ=FG.GradQ;
GradS=FG.GradS;
GradO=FG.GradO;

figure; 
%subplot(4,1,1)
subplot(221)
plot(Spos,field,'om-','LineWidth',1);set(gca,'Fontsize',12,'Fontweight','b')
xlabel('S [m]'); ylabel('B Field [T]');
grid on;

%subplot(4,1,2)
subplot(222)
plot(Spos,GradQ,'sk--','LineWidth',1);set(gca,'Fontsize',12,'Fontweight','b')
xlabel('S [m]'); ylabel('QP Gradient [T/m]');
grid on;

%subplot(4,1,3)
subplot(223)
plot(Spos,GradS,'*b:','LineWidth',1);set(gca,'Fontsize',12,'Fontweight','b')
xlabel('S [m]'); ylabel('SP Gradient [T/m²]');
grid on;

%subplot(4,1,4)
subplot(224)
plot(Spos,GradO,'.r-','LineWidth',1);set(gca,'Fontsize',12,'Fontweight','b')
xlabel('S [m]'); ylabel('OP Gradient [T/m³]');
grid on;

sgtitle(desc);




    
