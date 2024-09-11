function DAV = calcDA_smart_out( RING, nsteps, nturns, dp, z0, res, xmax, xmin, ymax,div)
%% Calculates Dynamic Aperture by tracking particles along radial lines on the XY plane until they are lost.
%% It does outward search
%% Mandatory input arguments
% RING: input lattice file; Tracking can be 6D or 4D as defined by the input lattice
% nsteps: Number of radial lines to search
% nturns: Numbers of turns to track
% dp: initial momentum or fixed momentum
% z0: initial longitudinal position [m]
% res: Resolution to find the DA
% xmax: limits of the region where to look for the DA [m]
% xmin: limits of the region where to look for the DA [m]
% ymax: limits of the region where to look for the DA [m]
%% Output parameters
% DAV: ((nsteps+1)x2) matrix containing DA border coordinates [m]
%% Usage example
% nsteps=24;nturns=2000;dp=0.0;z0=0.0;
% res=100e-6;
% xmax=0.012;xmin=-0.012;ymax=0.007
% DAV=calcDA_smart(RING,24,2000,0,0,100e-6,12e-3,-12e-3,7e-3)
%% History
% 2024/06/10 Saroj
angle_step=pi/nsteps;
DAV   = zeros(nsteps+1,2);
RRc=linspace(1e-3,11e-3,15);   % The radii where it starts searching
%% Generate the no. of turns to start with, as per the input:nturns
jj = 1;
while true
    den = div^(jj-1);
    numero = nturns / den;
    if numero >= 50 && numero <= 100
        jok = jj;
        break;
    end
    jj = jj + 1;
end
Nt = zeros(1, jok);
for jj = 1:jok
    Nt(jj) = round(nturns / div^(jj-1)); % default div = 2
end
Ntc = flip(Nt); % sequence of number of turns
parfor i=1:nsteps+1
    angle = (i-1)*angle_step;
    Nt=Ntc;
    RR=RRc;
    r0=RR(1);
    m=1;
    rstable=0;
    while m <= length(RR)
        x = r0*cos(angle);
        y = r0*sin(angle);
        if (x>xmax || x<xmin || y>ymax)
            loss=true;
        else
            [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',Nt(1));
        end
        if (loss)
            rstable=RR(m-1);
            break;
        else
            m=m+1;
            r0=RR(m);

        end
    end
    j=1;
    while j < length(Nt)
        x = rstable*cos(angle);
        y = rstable*sin(angle);
        if (x>xmax || x<xmin || y>ymax)
            loss=true;
        else
            [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',Nt(j+1));
        end
        if (loss)
            m=m-1;
            rstable=RR(m);
        else
            j=j+1;
            rstable=RR(m) ;
        end
    end
    %% Refinement in radius 
    r1=rstable;    % stable point
    r2=RR(m+1);    % Unstable point
    r_new=(r1+r2)/2;
    while abs(r2-r1) >= res
        x = r_new*cos(angle);
        y = r_new*sin(angle);
        if (x>xmax || x<xmin || y>ymax)
            loss=true;
        else
            [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',Nt(j));
        end
        if (loss)
            r2=r_new;
        else
            r1=r_new;
            rstable=r_new;
        end
        r_new=(r1+r2)/2;
    end
    DAV(i,:)=[rstable*cos(angle),rstable*sin(angle)];
    disp([' DAV(' num2str(i) ') = ' num2str(DAV(i,:)) ])
end