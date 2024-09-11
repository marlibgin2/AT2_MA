function DAV = calcDA_smart( RING, nsteps, nturns, dp, z0, res, xmax, xmin, ymax)
%% Calculates Dynamic Aperture by tracking particles alogm radial lines on the XY plane until they are lost.

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

RR=linspace(2e-3,11e-3,15);   % The radii where it starts searching

%% Generate the no. of turns to start with, as per the input:nturns 
jj = 1;
while true
    den = 2^(jj-1);
    numero = nturns / den;
    if numero > 50 && numero < 100
        jok = jj;
        break;
    end
    jj = jj + 1;
end
Nt = zeros(1, jok); 
for jj = 1:jok
    Nt(jj) = round(nturns / 2^(jj-1));
end
Nt = flip(Nt); % sequence of number of turns

%% Precompute maximum stable points for all radii for a small number of turns 
Rm   = zeros(nsteps+1,2);
%r0=RR(1);
parfor i=1:nsteps+1
    angle = (i-1)*angle_step;
    r0=RR(1);
    %Rm(i,:) = calcDA_radii( RING, r0,dp,z0,xmax, xmin, ymax);
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
     Rm(i,:)=[rstable m-1];
end 

%% parallel computation for all radii and maximum number of turns
parfor i=1:nsteps+1
    r_stable=0;
    angle = (i-1)*angle_step;
    j=1;
    Rmc=Rm;       % make a local copy of Rm to avoid warning
    r=Rmc(i,1);
    k=Rmc(i,2);
    RRc=RR;       % make a local copy of RR to avoid warning 
    while j < length(Nt)
        
        x = r*cos(angle);
        y = r*sin(angle);

        [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',Nt(j+1));

        if (loss)  
            k=k-1;
            r=RRc(k);   
        else
            j=j+1;
            r=RRc(k) ;
            r_stable=r;
        end
    end
   
    %% Refinement in radius after it searches for stable and unstable radius for maximum number of turns by binary search method

    r1=r_stable;    % stable point
    r2=RRc(k+1);    % Unstable point
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
            r_stable=r_new;
        end
        r_new=(r1+r2)/2
    end
    DAV(i,:)=[r_stable*cos(angle),r_stable*sin(angle)];
    
end

    function Rs = calcDA_radii(RING,r0,dp,z0,xmax, xmin, ymax)
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
        Rs=[rstable m-1];
    end
end  