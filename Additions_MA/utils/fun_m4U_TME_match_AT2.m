function penalty = fun_m4U_TME_match_AT2(X0)
global THERING

% X1 = [5.4 -1.1]

%TOTangle = mTME19_match(X0);
%TOTangle = mTME19_match_AT2(X0);

Rout = m4U_TME_match_AT2(X0);
THERING = Rout; 

% 
% (X0(1), X0(2), X0(3), X0(4),  X0(5),  X0(6),  ...
%       X0(7), X0(8), X0(9), X0(10), X0(11))
AAA = atsummary;
emix = AAA.naturalEmittance;
xi  = AAA.chromaticity;
nu   = AAA.tunes;
[~,lindata]=atlinopt6(Rout,1:length(Rout)+1);
dphi = lindata(29).mu(1:2)/2/pi; % AAA.tunes;

U0   = AAA.radiation*1e6; % energy loss in keV
tw   = twissring(THERING,0.001,1:length(THERING)+1);
beta = cat(1,tw.beta);
bx   = beta(:,1);
by   = beta(:,2);

disp(['emix = ' num2str(emix*1e12) ' (pm)'])
disp(['  U0 = ' num2str(U0) ' (keV)'])

disp(['  nu = ' num2str(nu) ''])
disp(['  xi = ' num2str(xi) ''])
disp([' cell dphi = ' num2str(dphi) ''])

stiff = sqrt(sum((xi./nu).^2));

disp(['  stiffness = ' num2str(stiff) ''])
%disp([' cell length = ' num2str(AAA.circumference) ' (m)'])


%%%%% reference function with control of chromaticity

%% fn. without control of chromaticity
%%% emix = emix * (emix/10e-12)^3;
bxi = bx(1);
byc = by(16);%by(15);
bxc = bx(16);
dphix = dphi(1);
dphiy = dphi(2);

%%%%%%penalty = (emix/33e-12)^4 + (bxi/13)^4 + (byc/13)^4 + (bxc/20)^4; %((TOTangle*280/pi-2)/.001)^4 ; 

penalty = (emix/75e-12)^4 + (bxi/9)^4 + (byc/8)^4 + (bxc/0.35)^4 + ...
          (xi(1)/20)^6 + (xi(2)/50)^6 + ...
          ((dphix - 3/7)/0.01)^6 + ((dphiy - 1/7)/0.01)^6; 
        %  (U0/20)^6 + (xi(1)/3.8)^6 + (20*TOTangle*180/pi-360)^6; %((TOTangle*280/pi-2)/.001)^4 ; 

if isnan(nu(1)*nu(2))
     penalty = 1e9;
end
if ~isreal(nu(1)) || ~isreal(nu(2)) 
     penalty = 1e9;
end
end%