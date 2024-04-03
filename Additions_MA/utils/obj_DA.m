function DA = obj_DA(X, rrr)
nlines = 7;
dpoffset = 0.0; 
npointsperline = 10;
nturns = 500;
xmax = 10e-3; ymax = 6e-3;

if 1 == 0
%rrr,nlines,npointsperline,nturns,xmax,ymax,dpoffset)
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


verbo = 0;
theta = linspace(-pi/2,pi/2,nlines);
for nn=1:nlines
    x(nn,:) = linspace(0,xmax*sin(theta(nn)),npointsperline);
    y(nn,:) = linspace(0,ymax*cos(theta(nn)),npointsperline);
end


if verbo > 0
    fprintf('tracking %d particles for %d turns at dp = %f%% ...\n',nlines*npointsperline,nturns,dpoffset*100)
end
    %[o4, fp4] = findorbit4(THERING,dpoffset,1:length(THERING)); % calculate closed orbit 

    
    try

    [o4]   = findorbit4(rrr,dpoffset,1:length(rrr)); % calculate closed orbit 
r      = zeros(6,nlines*npointsperline);
r(1,:) = x(:);
r(3,:) = y(:);
r(5,:) = dpoffset; ringpass(rrr,r,1);
r1     = ringpass(rrr,r,nturns,'KeepLattice'); %'reuse'
x1  = r1(1,:);
px1 = r1(2,:);
y1  = r1(3,:);
py1 = r1(4,:);

x2  = reshape(x1,nlines,npointsperline,nturns);
px2 = reshape(px1,nlines,npointsperline,nturns);
y2  = reshape(y1,nlines,npointsperline,nturns);
py2 = reshape(py1,nlines,npointsperline,nturns);

for nn = 1:nlines
    index = find(isnan(squeeze(x2(nn,:,end))));
    if ~isempty(index)
        if index(1)>1
            edge(nn) = index(1)-1;
        else
            edge(nn) = 1;
        end
    else
        edge(nn) = npointsperline;
    end
    x_da(nn) = x(nn,edge(nn));
    y_da(nn) = y(nn,edge(nn));
end

%[k,DA] = boundary(x_da',y_da');
wx=1.3; wy=0.7; 
DA     = -sum(sqrt((wx*x_da).^2+(wy*y_da).^2));
    catch
DA = -1e19; 
    end
disp(['DA = ' num2str(DA)])
end