function [penalty, Rout] = fun_medmax_match_AT2(X, Rin, graf, fitchro)
if nargin<2
    graf=0;
end
%global Rin

dph1    = [NaN NaN];
dph2    = [NaN NaN];
dph3    = [NaN NaN];
dph4    = [NaN NaN];

qfii    = findcells(Rin,'FamName','QFI');
qfoi    = findcells(Rin,'FamName','QFO');
qfmi   = findcells(Rin,'FamName','QFM');
qfendi = findcells(Rin,'FamName','QFEND');
qdendi = findcells(Rin,'FamName','QDEND');

sdi     = findcells(Rin,'FamName','SD');
sdendi  = findcells(Rin,'FamName','SDEND');
sfmi    = findcells(Rin,'FamName','SFM');
sfoi    = findcells(Rin,'FamName','SFO');
sfii    = findcells(Rin,'FamName','SFI');

dipi   = findcells(Rin,'FamName','DIP');
dipmi  = findcells(Rin,'FamName','DIPm');
for i=1:length(dipi)
    wdip(i) = Rin{dipi(i)}.PolynomB(2)/Rin{dipi(7)}.PolynomB(2);
end
for i=1:length(dipmi)
    wdipm(i) = Rin{dipmi(i)}.PolynomB(2)/Rin{dipmi(7)}.PolynomB(2);
end

VARi = {qfii; qfoi; qfmi; qfendi; qdendi};

for j = 1:5
for i=1:length(VARi{j})
    Rin{VARi{j}(i)}.PolynomB(2) = X(j);
    Rin{VARi{j}(i)}.K           = X(j);
end
end

clear VARi;
VARi = {dipi; dipmi};
W    = {wdip; wdipm};
for j = 1:2
    for i=1:length(VARi{j})
        Rin{VARi{j}(i)}.PolynomB(2) = X(j+5) * W{j}(i);
        Rin{VARi{j}(i)}.K           = X(j+5) * W{j}(i);
    end
end

VARi = {sdi; sdendi; sfmi; sfoi; sfii};
for j = 1:5
for i=1:length(VARi{j})
    Rin{VARi{j}(i)}.PolynomB(3) = X(j+7);
    Rin{VARi{j}(i)}.K           = X(j+7);
end
end
% verify physical quantities
%dpindex = findcells(THERING,'BendingAngle');

%
% force chromaticity to be 1/1
%
if fitchro==1
try
    disp('fitting chromaticity to [1,1]')
    Rin = atfitchrom(Rin,[1 1],'SFM','SDEND');
catch
    disp('cannot fit chromaticity ... mmachine unstable')
end
end
[M44, MS, orb] = findm44(Rin,0,1:length(Rin)+1);% [REFPTS,NE+1]);
%[M44, MS, orb] = findm44(rrr,0,dpindex);% [REFPTS,NE+1]);
cos_mu_x = (M44(1,1)+M44(2,2))/2;
cos_mu_y = (M44(3,3)+M44(4,4))/2;

sin_mu_x = sign(M44(1,2))*sqrt(-M44(1,2)*M44(2,1)-(M44(1,1)-M44(2,2))^2/4);
sin_mu_y = sign(M44(3,4))*sqrt(-M44(3,4)*M44(4,3)-(M44(3,3)-M44(4,4))^2/4);

ax = (M44(1,1)-M44(2,2))/2/sin_mu_x;
ay = (M44(3,3)-M44(4,4))/2/sin_mu_y;

bx = M44(1,2)/sin_mu_x;
by = M44(3,4)/sin_mu_y;

BX = squeeze((MS(1,1,:)*bx-MS(1,2,:)*ax).^2 + MS(1,2,:).^2)/bx;
BY = squeeze((MS(3,3,:)*by-MS(3,4,:)*ay).^2 + MS(3,4,:).^2)/by;

AX = -squeeze((MS(1,1,:)*bx-MS(1,2,:)*ax).*(MS(2,1,:)*bx-MS(2,2,:)*ax) + MS(1,2,:).*MS(2,2,:))/bx;
AY = -squeeze((MS(3,3,:)*by-MS(3,4,:)*ay).*(MS(4,3,:)*by-MS(4,4,:)*ay) + MS(3,4,:).*MS(4,4,:))/by;

if ~isreal(cos_mu_x) || ~isreal(cos_mu_y) || ~isreal(sin_mu_y) || ~isreal(sin_mu_y) || ...
        abs(cos_mu_x)>1 || abs(cos_mu_y)>1 || abs(sin_mu_y)>1 || abs(sin_mu_y)>1 || ...
        ~isreal(BX) || ~isreal(BY) || ~isreal(AX) || ~isreal(AY) || ...
        ~isreal(bx) || ~isreal(by) || ~isreal(ax) || ~isreal(ay) || ...
        ~isreal(MS)

    ex = 1e19;
    xi   = [1 1]*1e19;
    nu   = [1 1]*1e19;
    U0   = 1e19;
    bx   = ones(1,length(Rin)+1) * 1e19;
    by   = ones(1,length(Rin)+1) * 1e19;
    hx   = ones(1,length(Rin)+1) * 1e19;

else
    try
        AAA = ringpara(Rin); %atsummary;
        [ringdata,lindata]=atlinopt6(Rin,1:length(Rin)+1);
        AAA.tunes     = lindata(end).mu/2/pi;
        AAA.dph1      = lindata(66+30*1).mu/2/pi - lindata(66+30*0).mu/2/pi;
        AAA.dph2      = lindata(66+30*2).mu/2/pi - lindata(66+30*1).mu/2/pi;
        AAA.dph3      = lindata(66+30*3).mu/2/pi - lindata(66+30*2).mu/2/pi;
        AAA.dph4      = lindata(66+30*4).mu/2/pi - lindata(66+30*3).mu/2/pi;
    catch
        AAA.emittx = 1e19;
        AAA.chroms = [1 1]*1e19;
        AAA.tunes  = [1 1]*1e19;
        AAA.U0     = 1e19;
        AAA.dph1   = 1e19;
        AAA.dph2   = 1e19;
        AAA.dph3   = 1e19;
        AAA.dph4   = 1e19;
    end
    ex = AAA.emittx;
    xi   = AAA.chroms;
    nu   = AAA.tunes;
    U0   = AAA.U0; %radiation*1e6; % energy loss in keV
    dph1 = AAA.dph1;
    dph2 = AAA.dph2;
    dph3 = AAA.dph3;
    dph4 = AAA.dph4;

    try 
        tw   = twissring(Rin,0,1:length(Rin)+1,'chrom',0.0001);
        beta = cat(1,tw.beta);
        bx   = beta(:,1);
        by   = beta(:,2);
        eta  = cat(2,tw.Dispersion);
        hx   = eta(1,:);
    catch
        bx   = ones(1,length(Rin)+1) * 1e19;
        by   = ones(1,length(Rin)+1) * 1e19;
        hx   = ones(1,length(Rin)+1) * 1e19;    
    end
end



%%%%% reference function with control of chromaticity

%% fn. without control of chromaticity
%%% ex = ex * (ex/10e-12)^3;
bxi = bx(1); byi=by(1); bxmax = max(bx(45:207)); bymax=max(by(45:207));
hxi = hx(1);

bxU1 = bx(66+30*0); bxU2 = bx(66+30); bxU3 = bx(66+30*2); bxU4 = bx(66+30*3); bxU5 = bx(66+30*4); 
byU1 = by(66+30*0); byU2 = by(66+30); byU3 = by(66+30*2); byU4 = by(66+30*3); byU5 = by(66+30*4); 

bxS1 = bx(48+30*0); bxS2 = bx(48+30); bxS3 = bx(48+30*2); bxS4 = bx(144); bxS5 = bx(144+30*1); bxS6 = bx(144+30*2);

%          ((xi(1)-1)/0.5)^4 + ...
%          ((xi(2)-1)/0.5)^4 + ...
%           ((nu(1)-54.40)/6)^2 + ...
%           ((nu(2)-18.20)/6)^2 + ...

EX = 'exp((ex-150e-12)/1e-12)';
HX = '(hxi/1e-5)^4';
Bmax = '(bxmax/10)^2 + exp((bymax-20)/1)';
Bdif_i = '((bxi-byi)/3)^2';
By_i   = 'exp((byi-5)/0.2)';
NuRatio = '( ((nu(1)/nu(2))-3)/.3 )^4';
DPH1x    = '((dph1(1)-0.38)/0.02)^4';
DPH2x    = '((dph2(1)-0.36)/0.02)^4';
DPH3x    = '((dph3(1)-0.36)/0.02)^4';
DPH4x    = '((dph4(1)-0.38)/0.02)^4';
DPH1y    = '((dph1(2)-0.09)/0.01)^4';
DPH2y    = '((dph2(2)-0.11)/0.01)^4';
DPH3y    = '((dph3(2)-0.11)/0.01)^4';
DPH4y    = '((dph4(2)-0.09)/0.01)^4';

penalty =  eval(EX) +...  ((ex-131e-12)/10e-12)^4 +
           eval(HX) + ...  ;% + (bxi/9)^4;%(exp(ex-150e-12)/1e-12-1)+ ... %  
           eval(Bmax) + ... ; %  (bymax/20)^4;
           eval(Bdif_i) + ...
           eval(By_i) + ...
           eval(NuRatio) + ... 
           eval(DPH1x) + ...
           eval(DPH2x) + ...
           eval(DPH3x) + ...
           eval(DPH4x) + ...
           eval(DPH1y) + ...
           eval(DPH2y) + ...
           eval(DPH3y) + ...
           eval(DPH4y) ;
%           ((bxU1-bxU3)/0.2)^4 + ...
%           ((bxU2-bxU3)/0.2)^4 + ...
%           ((bxU4-bxU3)/0.2)^4 + ...
%           ((bxU5-bxU3)/0.2)^4 + ...
%           ((byU1-byU3)/0.2)^4 + ...
%           ((byU2-byU3)/0.2)^4 + ...
%           ((byU4-byU3)/0.2)^4 + ...
%           ((byU5-byU3)/0.2)^4 + ...
%          ((bxS2-bxS3)/2)^4 + ...
%          ((bxS5-bxS4)/2)^4 + ...
          

disp(['ex = ' num2str(ex*1e12) ' (pm) ... P = '  num2str( eval(EX) ) ]) 
disp(['  U0 = ' num2str(U0) ' (keV)'])

disp(['  nu = ' num2str(nu) '     ... P =  ' num2str(((nu(1)-54.40)/6)^2) ' ' num2str(((nu(2)-18.20)/6)^2) ]) 
disp(['  nu1/nu2 = ' num2str(nu(1)/nu(2)) '     ... P = ' num2str( eval(NuRatio)) ])
disp(['  xi = ' num2str(xi) ''])

stiff = sqrt(sum((xi./nu).^2));

disp(['  stiffness = ' num2str(stiff) ''])
disp(['  bxMAX = ' num2str(bxmax) ' byMAX = ' num2str(bymax) '   ... P = ' num2str(eval(Bmax))])
disp(['  bxi = ' num2str(bxi) ' byi = ' num2str(byi) '   ... PdiffBeta = ' num2str( ((bxi-byi)/3)^2 ) ]);
% disp(['bxU1 - bxU3  ... P = ' num2str(((bxU1-bxU3)/0.2)^4) ])
% disp(['bxU2 - bxU3  ... P = ' num2str(((bxU2-bxU3)/0.2)^4) ])
% disp(['bxU4 - bxU3  ... P = ' num2str(((bxU4-bxU3)/0.2)^4) ])
% disp(['bxU5 - bxU3  ... P = ' num2str(((bxU5-bxU3)/0.2)^4) ])
% disp(['byU1 - byU3  ... P = ' num2str(((byU1-byU3)/0.2)^4) ])
% disp(['byU2 - byU3  ... P = ' num2str(((byU2-byU3)/0.2)^4) ])
% disp(['byU4 - byU3  ... P = ' num2str(((byU4-byU3)/0.2)^4) ])
% disp(['byU5 - byU3  ... P = ' num2str(((byU5-byU3)/0.2)^4) ])
% disp(['bxS2 - bxS3  ... P = ' num2str(((bxS2-bxS3)/2)^4) ])
% disp(['bxS5 - bxS3  ... P = ' num2str(((bxS5-bxS3)/2)^4) ])
disp(['  hxi = ' num2str(hxi)  '        ... P = ' num2str((hxi/1e-5)^4) ])
disp(['  dph1 = ' num2str(dph1(1)) ' / ' num2str(dph1(2)) ' ... P = ' num2str(eval(DPH1x)) ' / '  num2str(eval(DPH1y)) ])
disp(['  dph2 = ' num2str(dph2(1)) ' / ' num2str(dph2(2)) ' ... P = ' num2str(eval(DPH2x)) ' / '  num2str( ((dph2(2)-0.10)/0.1)^4 )])
disp(['  dph3 = ' num2str(dph3(1)) ' / ' num2str(dph3(2)) ' ... P = ' num2str(eval(DPH3x)) ' / '  num2str(((dph3(2)-0.10)/0.1)^4 ) ])
disp(['  dph4 = ' num2str(dph4(1)) ' / ' num2str(dph4(2)) ' ... P = ' num2str(eval(DPH4x)) ' / '  num2str(eval(DPH4y)) ]); %num2str(((dph4(2)-0.10)/0.001)^4) ])         
         
% + ... %(byc/13)^4 + (bxc/20)^4 + ...
%         (U0/20)^6 + (xi(1)/3.8)^6; %((TOTangle*280/pi-2)/.001)^4 ; 

if isnan(nu(1)*nu(2))
     penalty = 1e9;
end
if ~isreal(nu(1)) || ~isreal(nu(2)) 
     penalty = 1e9;
end

if graf
    atplot(Rin); xlim([0 528/20])
end
Rout=Rin; 

end