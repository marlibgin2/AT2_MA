function [xx,zz]=atdynap(ring,nt,dpp,rfrac,varargin)
%ATDYNAP		Compute the dynamic aperture
%
%
%[XX,ZZ]=ATDYNAP(RING,NTURNS,DPP,RFRAC)
%[XX,ZZ]=ATDYNAP(...,'verbose',VALUE)
%
%XX,ZZ :	limit of the dynamic aperture (betatron amplitudes in m)
%RING :		Structure for tracking
%NTURNS:	Number of turns
%DPP :		Off-momentum value (default: 0)
%RFRAC :	Resolution of the grid for checking the stability
%			as a fraction of the maximum stable amplitude
%			(default: 0.02)
%
%'verbose': Provides text output if VALUE is true (default: true)   
[verboseFlag,~]=getoption(varargin,'verbose',true);

np=10;
rlist=0:0.001:0.1;
if nargin < 4, rfrac=0.02; end
if nargin < 3, dpp=0.0; end

if isnumeric(dpp)
clorb=[findorbit4(ring,dpp);dpp;0];
else
   clorb=findorbit6(ring);
end

t1=linspace(0,pi,2*np+3);
xpmax=ascan(ring,nt,clorb,0,rlist,verboseFlag);
zmax=ascan(ring,nt,clorb,0.5*pi,rlist,verboseFlag);
xmmax=ascan(ring,nt,clorb,pi,rlist,verboseFlag);
% 
% x1=[xpmax*ones(1,np+2) xmmax*ones(1,np+1)];
% z1=zmax*ones(1,2*np+3);
% tlist=atan2(sin(t1).*z1,cos(t1).*x1)';
% 
% rr=NaN(2*np+3,1);
% rr(1)=xpmax;
% rr(np+2)=zmax;
% rr(2*np+3)=xmmax;
% for i=[2:np+1 np+3:2*np+2]
%    rr(i)=ascan(ring,nt,clorb,tlist(i),rlist);
% end
% xx=rr.*cos(tlist);
% zz=rr.*sin(tlist);
slist=0.5:rfrac:2;
xx=NaN(2*np+3,1);
zz=xx;
for i=1:np+3
   [xx(i),zz(i)]=bscan(ring,nt,clorb,...
	  xpmax*cos(t1(i))*slist,zmax*sin(t1(i))*slist,verboseFlag);
end
for i=np+4:2*np+3
   [xx(i),zz(i)]=bscan(ring,nt,clorb,...
	  xmmax*cos(t1(i))*slist,zmax*sin(t1(i))*slist,verboseFlag);
end

function rmax=ascan(ring,nt,clorb,theta,rlist,verboseFlag)
for rr=rlist
   rin=clorb+[rr*cos(theta);0;rr*sin(theta);0;0;0];
   [dummy,lost]=ringpass(ring,rin,nt,'KeepLattice'); %#ok<ASGLU>
   if lost, break; end
   rmax=rr;
end
if verboseFlag
    fprintf('theta: %g, r: %g\n',theta,rmax);
end

function [xmax,zmax]=bscan(ring,nt,clorb,xlist,zlist,verboseFlag)
xmax = 0.0;
zmax = 0.0;
for i=1:length(xlist)
   rin=clorb+[xlist(i);0;zlist(i);0;0;0];
   [dummy,lost]=ringpass(ring,rin,nt,'KeepLattice'); %#ok<ASGLU>
   if lost, break; end
   xmax=xlist(i);
   zmax=zlist(i);
end
if verboseFlag
    fprintf('xm: %g, zm: %g\n',xmax,zmax);
end
