function [tau_p, tau_h, tau_v] = ibs_bane(RING,N,s_s,s_p,e_h,e_v)
% function [tau_p, tau_h, tau_v] = ibs_bane(RING,N,s_s,s_p,e_h,e_v)
%  Finds H-functions of RING and calculates IBS growth rates
%  using Karl Bane's approximated formulas (at every lattice point).
%
%  It calls FindOrbit4 and LinOpt which assume a lattice
%  with NO accelerating cavities and NO radiation
% 
% 2005, Christoph Steier

re = 2.8179409e-15;
c = 2.99792458e8;

global GLOBVAL;

L = length(RING);
spos = findspos(RING,1:L+1);
[lopt, tunes] = twissring(RING,0,1:L);

RINGLength = spos(L)+RING{L}.Length; 
betax = zeros(1,L);
betay = zeros(1,L);
alphax = zeros(1,L);
alphay = zeros(1,L);
gammax = zeros(1,L);
gammay = zeros(1,L);

%%% get twiss parameters
for i =1:L 
    betax(i)  = lopt(i).beta(1);
    betay(i)  = lopt(i).beta(2);
    alphax(i) = lopt(i).alpha(1);
    alphay(i) = lopt(i).alpha(2);
end
gammax = (1+alphax.*alphax)./betax;
gammay = (1+alphay.*alphay)./betay;    

betax(L+1)  = betax(1);
betay(L+1)  = betay(1);
alphax(L+1) = alphax(1);
alphay(L+1) = alphay(1);
gammax(L+1) = gammax(1);
gammay(L+1) = gammay(1);
tunes;


%%% get dispersion
orbit = findorbit4(RING,0.001,1:length(RING)+1)*1000;
etax = orbit(1,:);
etaxp= orbit(2,:);
etay = orbit(3,:);
etayp= orbit(4,:);

Hx = gammax.*etax.*etax + 2*alphax.*etax.*etaxp + betax.*etaxp.*etaxp;
Hy = gammay.*etay.*etay + 2*alphay.*etay.*etayp + betay.*etayp.*etayp;

%///////////////////////////////////////////////////////////////////////////////
% figure
% subplot(2,1,1)
% plot(spos,abs(betax),'b');
% hold on
% plot(spos,abs(betay),'r');
% A= axis;
% A(1) = 0;
% A(2) = RINGLength;
% axis(A);
% xlabel('s - position [m]');
% ylabel('\beta_x \beta_y [m]');
% grid on
% 
% title('beta-functions and H-functions');
% subplot(2,1,2)
% plot(spos,Hx,'b');
% hold on
% plot(spos,Hy,'r');
% A= axis;
% A(1) = 0;
% A(2) = RINGLength;
% axis(A);
% xlabel('s - position [m]');
% ylabel('H_x H_y [m]');
% grid on
% addlabel(datestr(now));
%/////////////////////////////////////////////////////////////////////////////////
gamm = (GLOBVAL.E0./511e3);

clog = log(gamm.^2.*sqrt(mean(abs(betay)).*e_v).*e_h./(re.*mean(abs(betax))));

sigma_H = (1./s_p.^2+Hx./e_h+Hy./e_v).^(-1/2);

a = sigma_H./gamm.*sqrt(betax./e_h);
b = sigma_H./gamm.*sqrt(betay./e_v);
alph=a./b;
g_bane=alph.^(0.021-0.044*log(alph)); %//// Hamed_ Bane2012

diffs=[diff(spos) 0];

% g_ave = 0;
% for loop = 1:length(diffs)
%     g_ave = g_ave+diffs(loop)*sigma_H(loop)*g_bane(a(loop)/b(loop)) ...
%         *(betax(loop)*betay(loop))^(-1/4);
% end
% g_ave = g_ave/spos(end);
%     
% tau_p = ((re.^2*c*N*clog)./(16.*gamm.^3.*e_h.^(3/4).*e_v.^(3/4).*s_s.*s_p.^3).*g_ave).^(-1);


tau_p = ((re.^2*c*N*clog)./(16.*gamm.^3.*e_h.^(3/4).*e_v.^(3/4).*s_s.*s_p.^3) ...
    .*sum(diffs.*sigma_H.*g_bane.*(betax.*betay).^(-1/4))./spos(end)).^(-1);  % *g_bane(a./b).*

tau_h = e_h./(s_p.^2.*sum(diffs.*Hx)./spos(end)).*tau_p;
tau_v = e_v./(s_p.^2.*sum(diffs.*Hy)./spos(end)).*tau_p;