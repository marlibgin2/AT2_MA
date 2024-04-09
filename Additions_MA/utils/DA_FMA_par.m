% -------------------------------
% calculate DA/NU diffusion rates 
% -------------------------------
function [x0pos, y0pos, nuxpos, nuypos, diffuvec] = DA_FMA_par(ring,filename,nturn,nx,ny,xmax,ymax)
pp = gcp; 

display_output = 1;
if nargin<2
    filename = 'fma';
    nturn=500;
    nx=20; 
    ny=20; 
    xmax = 7e-3;
    ymax = 7e-3;
end

d = '';                     % directory to save data
outfile = [filename '.out'];% filename

% -----------------------------------------
% define grid of initial input for tracking
% -----------------------------------------
r = zeros(6,nx*ny);
x = repmat(linspace(-xmax,xmax,nx),1,ny);
vec = repmat(linspace(1e-6,ymax,ny),nx,1);
y = vec(:).';
r(1,:) = x;
r(3,:) = y;
dx = (x(2)-x(1))/2;                  % x-step size
dy = (vec(1,2) - vec(1,1))/2;        % y-step size  

% --------------------------------------------------
% calculate DA diffusion rate for every grid element
% --------------------------------------------------
da_fma_fast_mach(ring, nturn, r, [d outfile],display_output);
[x0pos, y0pos, nuxpos, nuypos, diffuvec] = plot_fma_machine(nturn, [0 1 0 1], [d outfile], dx, dy,display_output);

delete(pp)

function da_fma_fast_mach(RING,  nturn, r, outfile, display_output)
vv = [];
syms vv;
v  = {};
syms v; 
jj = 0; syms jj;
    
nfreq = 1; % looking for 6 frequencies beside the tune: put 1 for plain FMA

% Initialize ringpass
ringpass(RING, [1e-4 0 0 0 0 0]');

% establish twiss/global parameters
dp = 1e-8;
[Twiss, tune, chrom] = twissring(RING, dp, 1:(length(RING)+1), 'chrom');

fid = fopen(outfile, 'w');
fprintf(fid, 'Dynamic Tracking with NAFF \n');
fclose(fid);

nturn2 = nturn/2;                % track the 1st half of turns ... 
fid = fopen(outfile, 'a');

% track a full "grid" of particles with ringpass
disp('Tracking particles....')
v={}; vv = zeros(6,nx*ny,nturn);
parfor i=1:nx*ny
    RR = r(:,i);
    v{i} = ringpass(RING, RR, nturn);
    %vv(1:6,i,1:nturn) = ringpass(RING, r(:,i), 1:nturn);
end
for i=1:nx*ny
    vv(1:6,i,1:nturn) = v{i}(1:6,1:nturn); 
end
clear v; v = vv; 

% v = ringpass(RING, r, nturn);
% v = reshape(v,6,size(r,2),nturn);
disp('Done')
nux = [];
nuy = [];
parfor i = 1:size(r,2)
    fvect = []; frequency=[]; amplitude=[]; phase=[]; 
    xf = squeeze(v(:,i,:));
    xfpeak=xf(1,nturn);
    zfpeak=xf(3,nturn);

    if(isnan(xfpeak) || isnan(zfpeak) || (xfpeak == 0) || (zfpeak == 0)) % with lost particles we save zeros ...
        if display_output
            fprintf('%s\n',' particle lost, printing all zeros');
        end
        nux(i) = 0; dnux(i) = 0;
        nuy(i) = 0; dnuy(i) = 0; 
        %fprintf(fid, '%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n', r(1,i),r(3,i), nux(i), nuy(i), dnux(i), dnuy(i));
    else   % if the particle has not been lost... we can save something meaningful
        if display_output
            disp('Survived all turns');
        end
        % A. First half of data
        mxhalf1 = mean(xf(:,1:nturn2)');
        xhalf = xf(:,1:nturn2)-mxhalf1'*ones(1,nturn2); % subtracting mean value
        %                xhalf = xf(:,1:nturn2);
       [freq amplitude phase] = calcnaff(xhalf(1,:), xhalf(2,:)); 
       [a b] = max(abs(amplitude)); 
       nf = length(freq); if (nf < nfreq), freq((nf+1):nfreq) = 0.0; end;
       nux(i) = abs(freq(b)/2/pi); 

       [freq amplitude phase] = calcnaff(xhalf(3,:), xhalf(4,:)); 
       [a b] = max(abs(amplitude)); 
       nf = length(freq); if (nf < nfreq), freq((nf+1):nfreq) = 0.0; end;
       nuy(i) = abs(freq(b)/2/pi); 
       
       %fprintf(fid, '%13.6E %13.6E %13.6E %13.6E ', r(1,i),r(3,i), nux(i) ,nuy(i));

       % B. Second half of data
       mxhalf2 = mean(xf(:,(nturn2+1):nturn)');
       xhalf = xf(:,(nturn2+1):nturn)-mxhalf2'*ones(1,nturn2); % subtracting mean value

       [freq amplitude phase] = calcnaff(xhalf(1,:), xhalf(2,:)); 
       [a b] = max(abs(amplitude)); 
       nf = length(freq); if (nf < nfreq), freq((nf+1):nfreq) = 0.0; end;
       nux_2 = abs(freq(b)/2/pi); 

       [freq amplitude phase] = calcnaff(xhalf(3,:), xhalf(4,:)); 
       [a b] = max(abs(amplitude)); 
       nf = length(freq); if (nf < nfreq), freq((nf+1):nfreq) = 0.0; end;
       nuy_2 = abs(freq(b)/2/pi); 
 
       dnux(i) = nux_2 - nux(i); dnuy(i) = nuy_2 - nuy(i); 
       %fprintf(fid, '%13.6E %13.6E\n', dnux(i), dnuy(i)); % nux_2 - nux, nuy_2 - nuy);
       if display_output
           fprintf('%13.6E %13.6E %13.6E %13.6E\n', r(1,i), r(3,i), nux(i), nuy(i));
       end
    end
end

for i = 1:length(nux)       
fprintf(fid, '%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n', r(1,i),r(3,i), nux(i), nuy(i), dnux(i), dnuy(i));
       %fprintf(fid, '%13.6E %13.6E %13.6E %13.6E ', r(1,i),r(3,i), nux(i) ,nuy(i));
       %fprintf(fid, '%13.6E %13.6E\n', dnux(i), dnuy(i)); % nux_2 - nux, nuy_2 - nuy);
end

fclose(fid);

end

end

function [x0pos, y0pos, nuxpos, nuypos, diffuvec] = plot_fma_machine(nturn, qrange, filestr, dx, dy,display_output)
%
% this MATLAB script plots the frequency map as computed by tracy-II
% in the file fmap.out
%
% nturn    = number of turns
% qrange   = [xmin xmax ymin ymax] for tunex vs tuney plot, including
%            integer tunes;
% filestr  = location of datafile and title string on the plots

% read in the data file
fmap=filestr;      % Tracy-II output

[x0 y0 fx fy dfx dfy] = textread(fmap,'%f %f %f %f %f %f','headerlines',1);

dnux=7e-3;%1e-3;
dnuy=7e-3;%1e-3;
dmx = dx;
dpx = dx;
dmy = dy;
dpy = dy;

% % reverting tune values to the interval [0, 1)
ntune=length(fx);
% for nt=1:ntune
%     if(fx(nt) <= 0)
%         fx(nt)=-fx(nt);
%     else
%         fx(nt)=1e0-fx(nt);
%     end
%     if(fy(nt) <= 0)
%         fy(nt)=-fy(nt);
%     else
%         fy(nt)=1e0-fy(nt);
%     end
% end

% building the arrays to be plotted
icounter=0;
for ny=1:ntune
    icounter=icounter+1;
    if (fx(icounter) ~= 0 || fy(icounter) ~= 0)
        % building a rectangular area in the x-y plot around each point with non zero tune
        x0grid(icounter)=x0(icounter);
        y0grid(icounter)=y0(icounter);
        
        % the rectangule is described in clockwise way from bottom left point
        x0pos(:,icounter) = ...
            [x0grid(icounter)-dmx; x0grid(icounter)-dmx; x0grid(icounter)+dpx; x0grid(icounter)+dpx];
        y0pos(:,icounter) = ...
            [y0grid(icounter)-dmy; y0grid(icounter)+dpy; y0grid(icounter)+dpy;
            y0grid(icounter)-dmy];
        % buiding rectangular area in the tune plot around each tune point
        nux(icounter)=fx(icounter);
        nuy(icounter)=fy(icounter);
        nuxpos(:,icounter) = ...
            [nux(icounter)-dnux; nux(icounter)-dnux; nux(icounter)+dnux; nux(icounter)+dnux];
        nuypos(:,icounter) = ...
            [nuy(icounter)-dnuy; nuy(icounter)+dnuy; nuy(icounter)+dnuy; nuy(icounter)-dnuy];
        % determining the diffusion coefficient
        delta_tune= dfx(icounter)^2+dfy(icounter)^2;
        if (delta_tune > 0)
            diffu(icounter) = sqrt(delta_tune);
            diffu(icounter) = log10(diffu(icounter)/(nturn/2));
            if (diffu(icounter) < (-10))
                diffu(icounter) = -10;
            end
        else
            diffu(icounter) = -10;
        end
        diffuvec(1:4,icounter) = diffu(icounter);
        diffuvec(1:4,icounter) = diffu(icounter);
    else
        x0pos(:,icounter)  = [100; 100; 100; 100];
        y0pos(:,icounter)  = [100; 100; 100; 100];
        nuxpos(:,icounter) = [0; 0; 0; 0];
        nuypos(:,icounter) = [0; 0; 0; 0];
        diffuvec(:,icounter) = [-10; -10; -10; -10];
        diffuvec(:,icounter) = [-10; -10; -10; -10];
    end
end

if display_output
    figure(35); clf
    fill(x0pos*1e3, y0pos*1e3, diffuvec);
    axis([-max(x0(:))*1e3, max(x0(:))*1e3, 0, max(y0(:))*1e3]);
    caxis([-10 -3]);
    grid on;
    hold on;
    shading flat;
    colormap('jet');
    colorbar;
    % tstr1 = ['nturn = ', num2str(nturn)];
    % title(tstr1);
    title('Dynamic Aperture at Injection Point')
    xlabel('Horizontal Amplitude (mm)');
    ylabel('Vertical Amplitude (mm)');
    
    figure(36); clf
    fill(nuxpos, nuypos, diffuvec);
    hold on;
    axis(qrange);
    caxis([-10 -3]);
    shading flat;
    colormap('jet');
    colorbar;
    
    plot_net(7,qrange(1),qrange(2),qrange(3),qrange(4));
    
    tstr2 = ['Frequency Map'];
    title(tstr2);
    xlabel('\nu_x');
    ylabel('\nu_y');
    
end
end