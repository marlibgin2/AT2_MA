function penalty = funTME19_match(X0)
global THERING

% ----------------------------------------------------------------------------------------------
% solutions from X3 to X11 are in file mTME10_optimisation.mat
% X12 = [ 4.0597   -3.5117    0.0035    0.0052    0.0949    0.0334  0.2606    0.0227    0.8259 ]
% ----------------------------------------------------------------------------------------------

TOTangle = mTME19_match(X0);

% 
% (X0(1), X0(2), X0(3), X0(4),  X0(5),  X0(6),  ...
%       X0(7), X0(8), X0(9), X0(10), X0(11))
AAA = atsummary;
emix = AAA.naturalEmittance;
xi  = AAA.chromaticity;
nu   = AAA.tunes;
U0   = AAA.radiation*1e6; % energy loss in keV
tw   = twissring(THERING,0.001,1:length(THERING)+1);
beta = cat(1,tw.beta);
bx   = beta(:,1);
by   = beta(:,2);

disp(['emix = ' num2str(emix*1e12) ' (pm)'])
disp(['  U0 = ' num2str(U0*20) ' (keV)'])

disp(['  nu = ' num2str(nu*20) ''])
disp(['  xi = ' num2str(xi*20) ''])

stiff = sqrt(sum((xi./nu).^2));

disp(['  stiffness = ' num2str(stiff) ''])
%disp([' cell length = ' num2str(AAA.circumference) ' (m)'])


%%%%% reference function with control of chromaticity

%% fn. without control of chromaticity
%%% emix = emix * (emix/10e-12)^3;
bxi = bx(1);
byc = by(73);%by(15);
bxc = bx(73);
%penalty = (emix/95e-12)^2 + (bxi/15)^4 + (byc/15)^4; 
%penalty = (emix/80e-12)^2 + (bxi/15)^4 + (byc/15)^4; 
%penalty = (emix/58e-12)^2 + (bxi/13)^4 + (byc/13)^4; 


%%%%%%penalty = (emix/33e-12)^4 + (bxi/13)^4 + (byc/13)^4 + (bxc/20)^4; %((TOTangle*280/pi-2)/.001)^4 ; 

penalty = (emix/30e-12)^4 + (bxi/13)^4 + (byc/13)^4 + (bxc/20)^4 + ...
          (U0/20)^6 + (xi(1)/3.8)^6 + (20*TOTangle*180/pi-360)^6; %((TOTangle*280/pi-2)/.001)^4 ; 

if isnan(nu(1)*nu(2))
     penalty = 1e9;
end
if ~isreal(nu(1)) || ~isreal(nu(2)) 
     penalty = 1e9;
end
end%