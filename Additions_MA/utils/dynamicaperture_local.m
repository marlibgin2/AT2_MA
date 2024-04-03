function [x0pos, y0pos, nuxpos, nuypos, diffuvec] = dynamicaperture_local(ring,filename,nturn,nx,ny,xmax,ymax)
%REFERENCE
%/gpfs/offline1/staff/common/marapo/DIAMONDII/LOCO_TESTS/MISALIGNMENT_afterMAC1/HGH/21-8-1
% rf = findcells(ring,'FamName','RF');,
% ring{rf}.Voltage = 0;

d = '';% directory to save data
outfile = [filename '.out'];% filename

display_output = 1;

nturn=132; %500;%1000;
nx=2; %40; %80;
ny=2; %40; %80;
xmax = 1e-4; %7e-3;
ymax = 1e-4; %7e-3;
r = zeros(6,nx*ny);
x = repmat(linspace(-xmax,xmax,nx),1,ny);
vec = repmat(linspace(1e-6,ymax,ny),nx,1);
y = vec(:).';
r(1,:) = x;
r(3,:) = y;
dx = (x(2)-x(1))/2;
dy = (vec(1,2) - vec(1,1))/2;

da_fma_fast_mach(ring, nturn, r, [d outfile],display_output);
% plot_fma_machine(nturn, [0.15 0.25 0.28 0.38], [d outfile], dx, dy);
[x0pos, y0pos, nuxpos, nuypos, diffuvec] = plot_fma_machine(nturn, [0 1 0 1], [d outfile], dx, dy,display_output);



function da_fma_fast_mach(THERING,  nturn, r, outfile,display_output)
% Usage:
% Tracking for dynamic aperture
%
% Input:
%  nloc       -- starting location for tracking. If 0, beginning of lattice
%  nturn      -- number of tracking turns
%  r          -- grid of particles to track
%  outfile    -- name of the o
% utput file in "tstring";

n = 4;
NAFF_MIN = nturn/2;
nwin = 1;
nfreq = 1; % looking for 6 frequencies beside the tune: put 1 for plain FMA
%tuneguessx =-0.225; % guess of Qx
%tuneguessy =-0.363; % guess of Qy
tuneguessx =-0.1; % guess of Qx
tuneguessy =-0.2; % guess of Qy
eps=0.02;           % allowed distance from the tune

% Initialize ringpass
ringpass(THERING, [1e-4 0 0 0 0 0]');

dp = 1e-8;
[Twiss, tune, chrom] = twissring(THERING, dp, 1:(length(THERING)+1), 'chrom');

fid = fopen(outfile, 'w');
fprintf(fid, 'Dynamic Tracking with NAFF: AT-Tracy 2002\n');
fclose(fid);

nturn2 = nturn/2;
fid = fopen(outfile, 'a');

% do tracking as a block
%tic
disp('Tracking particles....')
v = ringpass(THERING, r, nturn);
v = reshape(v,6,size(r,2),nturn);
disp('Done')
%toc

%tic
disp('Calculating FMAP and DA...')
% Make momentum and vertical as outer loops
for i = 1:size(r,2)
    xf = squeeze(v(:,i,:));
    % getting a NaN if the collimator is hit ***
    xfpeak=xf(1,nturn);
    zfpeak=xf(3,nturn);
    % NAFF (notice that a column vector is needed as input)
    if(isnan(xfpeak) | isnan(zfpeak) | (xfpeak == 0) | (zfpeak == 0))
        if display_output
            fprintf('%s\n',' particle lost, printing all zeros');
        end
        fprintf(fid, '%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n', r(1,i),r(3,i), 0, 0, 0, 0);
    else
        if display_output
            disp('Survived all turns');
        end
        % A. First half of data
        mxhalf1 = mean(xf(:,1:nturn2)');
        xhalf = xf(:,1:nturn2)-mxhalf1'*ones(1,nturn2); % subtracting mean value
        %                xhalf = xf(:,1:nturn2);
% % %         clear calcnaff;
        fvect = naff2(xhalf(1,:), [0.05 0.48]);
        itmp = 0;

        nf = length(fvect); if (nf < nfreq), fvect((nf+1):nfreq) = 0.0; end;
  
        nux = fvect(1);
% % %         clear calcnaff;
        fvect = naff2(xhalf(3,:), [0.05 0.48]);
        itmp = 0;
        % loop fino a itmp < 5 per trovare le frequenze nel piano
        % verticale
        nf = length(fvect); if (nf < nfreq), fvect((nf+1):nfreq) = 0.0; end;
        
        nuy = fvect(1);
        fprintf(fid, '%13.6E %13.6E %13.6E %13.6E ', r(1,i),r(3,i), nux ,nuy);
        
        % B. Second half of data
        mxhalf2 = mean(xf(:,(nturn2+1):nturn)');
        xhalf = xf(:,(nturn2+1):nturn)-mxhalf2'*ones(1,nturn2); % subtracting mean value
        %                xhalf = xf(:,(nturn2+1):nturn);
% % %         clear calcnaff;
        fvect = naff2(xhalf(1,:), [0.05 0.48]);
        itmp = 0;
        % loop fino a itmp < 5 per trovare le frequenze nel piano
        % orizzontale

        nf = length(fvect); if (nf < nfreq), fvect((nf+1):nfreq) = 0.0; end;
        
        nux_2 = fvect(1);
% % %         clear calcnaff;
        fvect = naff2(xhalf(3,:), [0.05 0.48]);
        itmp = 0;

        nf = length(fvect); if (nf < nfreq), fvect((nf+1):nfreq) = 0.0; end;
        
        nuy_2 = fvect(1);
        fprintf(fid, '%13.6E %13.6E\n', nux_2 - nux, nuy_2 - nuy);
        if display_output
            fprintf('%13.6E %13.6E %13.6E %13.6E\n', r(1,i), r(3,i), nux, nuy);
        end
    end
end; 

fclose(fid);
%toc


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
    figure(25); clf
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
    
    figure(26); clf
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