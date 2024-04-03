function [DA_Area, DA]= obj_modelDA(X, Rin, NT, fitchro)
if nargin < 3
    NT=100;
end
if 1  == 1
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
end

clear VARi; 

oxxoi = findcells(Rin,'FamName','oxxo');
oxyoi = findcells(Rin,'FamName','oxyo');
oyyoi = findcells(Rin,'FamName','oyyo');
VARi = {oxxoi; oxyoi; oyyoi};
% zero all the octupoles )in force from 24-11-2023)
noOCT=1;
if noOCT==1
for j=1:3
    for i=1:length(VARi{j})
        Rin{VARi{j}(i)}.PolynomB(4) = 0;
    end
end
end
%
% force chromaticity to be 1/1 (as it happens in the optimisation)
%
try
    disp('fitting chromaticity to [1,1]')
    Rin = atfitchrom(Rin,[1 1],'SFM','SDEND'); 
catch
    disp('cannot fit chromaticity ... mmachine unstable')
end

dp = 1e-8;
[Twiss, tune, chrom] = twissring(Rin, dp, 1:(length(Rin)+1), 'chrom');

r0      = 0.0125; % 0.02
nsteps  = 14; %9
nturns  = NT;
dp      = 0.0;
res     = 0.00035/4; % 0.0005

% modelDA( r0, nsteps, nturns, dp, res)
% Evalutes the Dynamic Aperture by tracking
% required arguments
% r0:     Inital amplitude X to search ~ 0.015 m
% nsteps: Number of points to find     ~ 20
% nturns: Numbers of turns to track    ~ 256
% dp:     Energy deviation             ~ 0.0%
% res:    Resolution to find the DA    ~ 0.0005 m
% Returns the Dynamic aperture in DA and the Data for the trackingrrr


% Written by M. Munoz

r_stable=0;
angle_step=pi/nsteps;

angle=0;
%global THERING;
r=r0;
try
    %Evaluate the Chromatic orbit
    twiss=  gettwiss(Rin, 0.0);
    x0=twiss.etax(1)*dp;
    %Check that the chromatic orbit is stable
    [T, loss]=ringpass(Rin,[x0 0.0 0 0.0 dp 0.0]',nturns);
    repeat=0;
    while repeat <2
    if (loss)
        disp('The chromatic closed orbit is not stable. No DA found');
        repeat = repeat+1; 
        DA(1,1)=0;
        DA(1,2)=0;
        Data=0;
    else
        repeat = 1000; 
        for i=1:nsteps+1
            look=true; r_stable=0;
            %%% fprintf('Tracing step %ld of %ld\n', i, nsteps);
            cnt=0;
            while (look)
                cnt=cnt+1;

                x= x0+r*cos(angle);
                y= r*sin(angle);
                %%%if (mod(cnt,1)==0)
                    %%%disp(['cnt = ' num2str(cnt) ' x= '  num2str(x) ' y= ' num2str(y) ])
                %%%end
                if cnt>=30
                    look=false;
                end
                [T, loss]=ringpass(Rin,[x 0.0 y 0.0 dp 0.0]',nturns);
                %%% fprintf('%s %d %d   \n','Tracked',r, angle);

                if (loss)
                    if ((r-r_stable) < res)
                        look=false; % we have found the limit of the DA
                        DA(i,1)=r_stable*cos(angle);
                        DA(i,2)=r_stable*sin(angle);
                        r=r_stable;
                        disp(['>>>> xf= '  num2str(DA(i,1)) ' yf= ' num2str(DA(i,2)) ])
                    else
                        r=(r+r_stable)/2; %2;
                    end
                else
                    r_stable=r;
                    r=r*1.1;%2
                    Data{i}=T;
                end
            end
            angle=angle+angle_step;
        end
    end
    end
    %SumDA = -sum(sqrt((DA(:,1).^2+(DA(:,2)).^2)));
    DA_Area = -calcDA_Area(DA);

catch
    DA_Area = 1e19;

end

end

function DA_Area = calcDA_Area(DA)
a = [];
A = [0, 0];
for i = 1:length(DA)-1
    B = DA(i,  :);
    C = DA(i+1,:);
    a(i) = (A(1) * (B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))/2;
end
DA_Area = sum(a);

end