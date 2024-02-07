% -------------------------------
% calculate DA/NU diffusion rates 
% using parallel capabilites
%
% from an original file used at 
% Diamond Light Source
%
% makes use of at-code calcnaff
% 
% adapted to work on clusters
% (Aurora, local... ) at MAXIV
%
% MA 01012024
% -------------------------------

%function [x0pos, y0pos, nuxpos, nuypos, diffuvec, WA] = DA_NU_FMA_par(ring,filename,nturn,nx,ny,xmax,ymax)
function [x0pos, y0pos, nuxpos, nuypos, diffuvec, WA] = DA_NU_FMA_par(ring,filename,nturn,nx,ny,xmax,ymax,varargin)

% -------------
% default input
% -------------
if nargin<3
    nturn    = 2048; %1056; 
    nx       = 50;
    ny       = 32;
    xmax     = 12e-3;
    ymax     = 7e-3;
end  
if nargin<2
    filename = 'fma';
end
display_output = 0; % verbosity flag  
graf           = 0;
calc           = 0; % new fresh calculation
controlgcp     = 0; % use cluster
filo           = 0; % decide whether to produce an output file 
CLUname        ='';
CLUnodes       = 0;

for ik = length(varargin):-1:1
    if strcmpi(varargin{ik},'verbosity')
        display_output = 1; % verbosity flag... 
        varargin(ik)   = [];        
    elseif strcmpi(varargin{ik},'newcalc')
        calc           = 1; % 1= new fresh calculation
        varargin(ik)   = [];        
%     elseif strcmpi(varargin{ik},'use_cluster')
%         controlgcp     = 1; % 1=use cluster
%         varargin(ik)   = [];     
    elseif strcmpi(varargin{ik},'aurora cluster')
        controlgcp     = 1; % 1=use cluster
        CLUname        = 'aurora R2022a'; % use maxiv Aurora cluster
        varargin(ik)   = [];        
    elseif strcmpi(varargin{ik},'local cluster')
        controlgcp     = 1; % 1=use cluster
        CLUname        = 'local'; % use maxiv Aurora cluster
        varargin(ik)   = [];        
    elseif strcmpi(varargin{ik},'cluster nodes')
        CLUnodes       = varargin{ik+1}; % use maxiv Aurora cluster
        varargin(ik+1)   = []; 
        varargin(ik) = [];         
    elseif strcmpi(varargin{ik},'graphic_plots')
        graf           = 1; % 1=use cluster
        varargin(ik)   = [];        
    elseif strcmpi(varargin{ik},'file_output')
        filo           = 1; % 1=use cluster
        varargin(ik)   = [];        
    end
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
% calc       = 1; % 1= new fresh calculation
% controlgcp = 1; % 1=use cluster
% graf       = 1; % 1=final plots / 0=no final plots
% filo       = 0; % decide whether to produce an output file 
if graf == 1
    filo = 1; % for graphic plots filo is needed
end

if calc==1
    if controlgcp==1
        c = parcluster;
        % next line is optional at MAX IV
        c.AdditionalProperties.AccountName = 'any-virtual-account-name';
        % 6 hour walltime
        c.AdditionalProperties.WallTime = '06:00:00';
        % hyperthreading enabled
        c.NumThreads = 6;
        c.saveProfile;
        %%%%parpool('aurora R2022a',56) %%% IMPORTANT REFERENCE ---
        %%%%parpool('aurora R2022a',56) %%% TEST with reduced nodes 
        %%%%parpool('local',12)
        if strcmpi(CLUname,'local')
            CLUnodes=ge(CLUnodes,12)*12+lt(CLUnodes,12)*CLUnodes; % local cluster max nodes is 12! 
        end
        parpool(CLUname, CLUnodes)
        pp = gcp; 
    end
    WA = da_fma_fast_mach(ring, nturn, r, [d outfile], display_output);    
end
if controlgcp==1
    delete(pp)
end

x0pos=[]; y0pos=[]; nuxpos=[]; nuypos=[]; diffuvec=[]; 
if graf==1
    [x0pos, y0pos, nuxpos, nuypos, diffuvec] = plot_fma_machine(nturn, [0 1 0 1], [d outfile], dx, dy, graf);
    %TuneEnergyDependence(ring); Ae=[];
    disp('calculating the tune-energy dependence ...')
    [qx, qy] = TuneEnergyDependence_tracking(ring,nturn,17,-4e-2,+4e-2);
    figure(36); hold on; axis([0 0.5 0 0.5]);
    plot(qx(10:17),qy(10:17),'ro','MarkerFaceColor','r')
    plot(qx(9),  qy(9),  'ko','MarkerFaceColor','k')
    plot(qx(1:8),  qy(1:8),  'co','MarkerFaceColor','c')
end


function WA = da_fma_fast_mach(RING,  nturn, r, outfile, display_output)
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
nturn2 = nturn/2;     % track the 1st half of turns ... 

if filo == 1
fid = fopen(outfile, 'w');
fprintf(fid, 'Dynamic Tracking with NAFF \n');
fclose(fid);
fid = fopen(outfile, 'a');
end
           


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
dnuxt=dnux; dnuyt=dnuy;
dnuxt(isnan(dnuxt))=0;
dnuyt(isnan(dnuyt))=0;
Sdnux = sum(abs(log10(dnuxt(dnuxt>0))));
Sdnuy = sum(abs(log10(dnuyt(dnuyt>0))));


% WA is the weighted normalized area of the diffusion plot
% where 10 is the ceiling of the nux,y diffusion
% consider the DA graph nx*ny pixels. If the whole plot were
% filled with "no-diffusion" pixels (dark blu dots) the 
% area would be nx*ny*10
% WA is therefeore the overall (sum) diffusion divided by (nx*ny*10)
WA    = (Sdnux + Sdnuy)/nx/ny/10; % NB: 10 is the ceiling of the nux,y diffusion !

if filo ==1
    for i = 1:length(nux)
        fprintf(fid, '%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n', r(1,i),r(3,i), nux(i), nuy(i), dnux(i), dnuy(i));
        %fprintf(fid, '%13.6E %13.6E %13.6E %13.6E ', r(1,i),r(3,i), nux(i) ,nuy(i));
        %fprintf(fid, '%13.6E %13.6E\n', dnux(i), dnuy(i)); % nux_2 - nux, nuy_2 - nuy);
    end
    fclose(fid);
end
end

end

function [x0pos, y0pos, nuxpos, nuypos, diffuvec] = plot_fma_machine(nturn, qrange, filestr, dx, dy, graf)
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

dnux=3.5e-3; %7e-3;%1e-3;
dnuy=3.5e-3; %7e-3;%1e-3;
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
        
        % the rectangle is described in clockwise way from bottom left point
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

if graf
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
    title('Dynamic Aperture at mid-STR01')
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
    
    res_ord=5; %7
    plot_net(res_ord,qrange(1),qrange(2),qrange(3),qrange(4));
    
    tstr2 = ['Tune Footprint'];
    title(tstr2);
    xlabel('\nu_x');
    ylabel('\nu_y');
    
end
end