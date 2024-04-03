function OBJ = OBJfcn2(X, rrr)
%global rrr

% -----------------------------------------
% alter the lattice according to parameters
% -----------------------------------------
rrr  = alter_lattice(X, rrr);


%
% force chromaticity to be 1/1
%
changechro=1;
if changechro == 1
    try
        disp('OBJ fitting chromaticity to [1,1]')
        rtemp = atfitchrom(rrr,[1 1],'SFM','SDEND');
         
        rrr=rtemp;  
    catch
        disp('cannot fit chromaticity ... mmachine unstable')
    end
end

try
    AAA = ringpara(rrr); %atsummary;
    if isnan(AAA.emittx)
        error('no emittance -- retry')
    end
catch
    AAA.emittx = 1e19;
end
emix = AAA.emittx;

if isnan(emix)
    emix = 1e19;
end


r0      = 0.02;
nsteps  = 9;
nturns  = 500;
dp      = 0.0;
res     = 0.0005;

r_stable=0;
angle_step=pi/nsteps;

angle=0;
%global THERING;
r=r0;
try
    %Evaluate the Chromatic orbit
    twiss=  gettwiss(rrr, 0.0);
    x0=twiss.etax(1)*dp;
    %Check that the chromatic orbit is stable
    [T, loss]=ringpass(rrr,[x0 0.0 0 0.0 dp 0.0]',nturns);
    if (loss)
        disp('The chromatic closed orbit is not stable. No DA found');
        DA(1,1)=0;
        DA(1,2)=0;
    else
        for i=1:nsteps+1
            look=true; r_stable=0;
            %%% fprintf('Tracing step %ld of %ld\n', i, nsteps);
            cnt=0;
            while (look)
                cnt=cnt+1;

                x= x0+r*cos(angle);
                y= r*sin(angle);
                % % %             if (mod(cnt,1)==0)
                % % %                 disp(['cnt = ' num2str(cnt) ' x= '  num2str(x) ' y= ' num2str(y) ])
                % % %             end
                if cnt>=30
                    look=false;
                    if r == 0
                        r = 1e-12;
                    end
                end
                [T, loss]=ringpass(rrr,[x 0.0 y 0.0 dp 0.0]',nturns);
                %%% fprintf('%s %d %d   \n','Tracked',r, angle);

                if (loss)
                    if ((r-r_stable) < res)
                        look=false; % we have found the limit of the DA
                        DA(i,1)=r_stable*cos(angle);
                        DA(i,2)=r_stable*sin(angle);
                        r=r_stable;
                    else
                        r=(r+r_stable)/2;
                    end
                else
                    r_stable=r;
                    r=r*2;
                    Data{i}=T;
                end
            end
            angle=angle+angle_step;
        end
    end
    DA_Area = -calcDA_Area(DA);

catch
    DA_Area = 1e19;

end

disp(['emix = ' num2str(emix) ' DA_Area = ' num2str(DA_Area)])
OBJ = [emix, DA_Area];

end

function DA_Area = calcDA_Area(DA)
a = [];
A = [0, 0];
for i = 1:length(DA)-1
    B = DA(i,  :);
    C = DA(i+1,:);
    a(i) = (A(1) * (B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))/2;
end
W = 1 - abs(DA(1,1) + DA(end,1))/(DA(1,1) - DA(end,1)); % weight to balance the X-asymmetry observed in some cases
DA_Area = sum(a) * W;
end

function rrr = alter_lattice(X, rrr)

if ~isempty(X)

    qfii    = findcells(rrr,'FamName','QFI');
    qfoi    = findcells(rrr,'FamName','QFO');
    qfmi   = findcells(rrr,'FamName','QFM');
    qfendi = findcells(rrr,'FamName','QFEND');
    qdendi = findcells(rrr,'FamName','QDEND');

    sdi     = findcells(rrr,'FamName','SD');
    sdendi  = findcells(rrr,'FamName','SDEND');
    sfmi    = findcells(rrr,'FamName','SFM');
    sfoi    = findcells(rrr,'FamName','SFO');
    sfii    = findcells(rrr,'FamName','SFI');

    dipi   = findcells(rrr,'FamName','DIP');
    dipmi  = findcells(rrr,'FamName','DIPm');
    for i=1:length(dipi)
        wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
    end
    for i=1:length(dipmi)
        wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
    end

    VARi = {qfii; qfoi; qfmi; qfendi; qdendi};

    for j = 1:5
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(2) = X(j);
            rrr{VARi{j}(i)}.K           = X(j);
        end
    end

    clear VARi;
    VARi = {dipi; dipmi};
    W    = {wdip; wdipm};
    for j = 1:2
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(2) = X(j+5) * W{j}(i);
            rrr{VARi{j}(i)}.K           = X(j+5) * W{j}(i);
        end
    end

    VARi = {sdi; sdendi; sfmi; sfoi; sfii};
    for j = 1:5
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(3) = X(j+7);
            rrr{VARi{j}(i)}.K           = X(j+7);
        end
    end

end
end