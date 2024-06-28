function TLT = calcTLT_raw(varargin)
% Calculates Touscheck lifetime following Piwinskys formalism
% Based on the AT2.0 TouschekPiwinskiLifeTime function, removing all
% optics calculations for speed
% 
%  momentum aperture column vector,  Nx2 column array with positive and
%                                 negative Local Momentum Apertures
%  latdat: linear lattice data structure                                 
%  current per bunch in A,        scalar
%  emittancex,                    default: atx modemittance(1)   scalar
%  emittancey,                    default: emittancex/2		     scalar
%  integration_method,            default: 'integral', may be: 'integral', 'quad', 'trapz')
%  energy_spread,                 scalar
%  bunch_length,	              scalar
%  AbsTol
%  RelTol
%			  
% OUTPUT
%
%  Tl  Lifetime in seconds 1/Tl=sum(contributionsTL.*L)/sum(L);
%

e0 = PhysConstant.elementary_charge.value; %1.60217646e-19; %Coulomb
r0 = PhysConstant.classical_electron_radius.value; %2.817940327e-15; %m %  classical electron radius
spl = PhysConstant.speed_of_light_in_vacuum.value; %299792458; % speed of ligth

%% Input argument parsing
[lindata,lmap] = getargs(varargin,[],[]);
Ib = getoption(varargin,'Ib',0.5/176);
verboselevel  = getoption(varargin,'verbose',0);
circumference = getoption(varargin,'circumference',528);
E0 = getoption(varargin,'energy',3.0e9);
emitx = getoption(varargin,'emitx', 328E-12/(1+0.025));
emity = getoption(varargin,'emity', 328E-12*0.025/(1+0.025));
sigp  = getoption(varargin,'sigp',7.8E-4);
sigs  = getoption(varargin,'sigs',0.010);
integrationmethod = getoption(varargin,'integrationmethod','integral');
abstol = getoption(varargin, 'AbsTol', 1.0e-16);
reltol = getoption(varargin, 'Relol', 1.0e-16);
tol = {'AbsTol', abstol, 'RelTol', reltol};

if (verboselevel>0)
    fprintf('*** Calculating Touschek lifetime with: \n');
    fprintf('emitx: %.3e [m]\n', emitx);
    fprintf('emity: %.3e [m]\n', emity);
    fprintf('energy spread: %.3e\n', sigp);
    fprintf('bunch length:  %.5g [m]\n', sigs);
    fprintf('integration method: "%s"\n', integrationmethod);
end
%% Preamble and calculations
% if dpp is a scalar assume constant momentum aperture.
%if numel(dpp)==1
%    dpp=dpp*ones(size(positions'));
%end
nps=size(lmap,1);
dppinput=lmap;
Tlcol=zeros(1,size(dppinput,2));
Nb = Ib/(spl/circumference)/e0; %Number of particle per bunch.
mass_el_ev=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;
relgamma = E0/mass_el_ev;%0.510998928e6;
relbeta=sqrt(1-1./relgamma.^2);

aaa=cat(1,lindata.alpha);
bbb=cat(1,lindata.beta);
ddd=cat(2,lindata.Dispersion)';

bx=bbb(:,1); % betax
by=bbb(:,2); % betay
Dx=ddd(:,1);
Dy=ddd(:,3);

ax=aaa(:,1);
ay=aaa(:,2);
Dpx=ddd(:,2);
Dpy=ddd(:,4);

sigxb=sqrt(emitx.*bx);
sigyb=sqrt(emity.*by);

sigx=sqrt(emitx.*bx+sigp.^2.*Dx.^2);
sigy=sqrt(emity.*by+sigp.^2.*Dy.^2);

Dtx=Dx.*ax+Dpx.*bx;%  % alpha=-b'/2
Dty=Dy.*ay+Dpy.*by;%

sigp2=sigp.^2;
Dx2=Dx.^2;
Dy2=Dy.^2;
Dtx2=Dtx.^2;
Dty2=Dty.^2;
sigxb2=sigxb.^2;
sigyb2=sigyb.^2;

sighinv2=1./(sigp2) +(Dx2+Dtx2)./(sigxb2) + (Dy2+Dty2)./(sigyb2);
sigh=sqrt(1./sighinv2);

B1=1./(2.*(relbeta.^2).*(relgamma.^2)).*( (bx.^2./(sigxb.^2)).*(1-(sigh.^2.*Dtx.^2./(sigxb.^2))) + (by.^2./(sigyb.^2)).*(1-(sigh.^2.*Dty.^2./(sigyb.^2))));
B2sq=1./(4.*(relbeta.^4).*(relgamma.^4)).*((bx.^2./(sigxb.^2)).*(1-(sigh.^2.*Dtx.^2./(sigxb.^2)))-(by.^2./(sigyb.^2)).*(1-(sigh.^2.*Dty.^2./(sigyb.^2)))).^2+(sigh.^4.*bx.^2.*by.^2.*Dtx.^2.*Dty.^2)./((relbeta.^4).*(relgamma.^4).*sigxb.^4.*sigyb.^4);
B2=sqrt(B2sq);
contributionsTL=zeros(size(dppinput));

for dppcolnum=1:size(dppinput,2)
    dpp=dppinput(:,dppcolnum);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%
    %%%%%%%% From here calculation takes place.
    %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    um=relbeta.^2*dpp.^2;
    %  em=bx.^2.*sigx.^2./(relbeta.^2.*relgamma.^2.*sigxb.^2.*sigtx2).*um;
    
    val=zeros(size(B1));
    km=atan(sqrt(um));
    FpiWfact=sqrt(pi.*(B1.^2-B2.^2)).*um;
    
    for ii=1:nps
        % choose integration method
        switch integrationmethod
            case 'infiniteintegral'
                val(ii)= integral(@(u)TLT_IntPiw(u,um(ii),B1(ii),B2(ii)),um(ii),Inf, tol{:}); %,...um(ii)*1e4
                
            case 'integral'
                val(ii) = integral(@(k)TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii)),km(ii),pi/2, tol{:}); %,...,'Waypoints',evalpoints um(ii)*1e4
                
            case 'quad'
                val(ii)= quad(@(k)TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii)),km(ii),pi/2); %,...,'Waypoints',evalpoints um(ii)*1e4
                
            case 'trapz'
                k=linspace(km(ii),pi/2,10000);
                val(ii)= trapz(k,TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii))); %,...,'Waypoints',evalpoints um(ii)*1e4
                
            otherwise % use default method quad (backward compatible)
                val(ii)=integral(@(k)TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii)),km(ii),pi/2, tol{:}); %,...,'Waypoints',evalpoints um(ii)*1e4
                
        end
    end
    
    
    
    switch integrationmethod
        case 'infiniteintegral'
            frontfact=(r0.^2.*spl.*Nb)./(8.*pi.*(relgamma.^2).*sigs.*sqrt(...
                (sigx.^2).*(sigy.^2)-sigp.^4.*Dx.^2.*Dy.^2).*um).*FpiWfact;
            
        case {'integral' 'quad' }
            frontfact=(r0.^2.*spl.*Nb)./(8.*pi.*(relgamma.^2).*sigs.*sqrt(...
                (sigx.^2).*(sigy.^2)-sigp.^4.*Dx.^2.*Dy.^2).*um).*2.*FpiWfact;
        case {'trapz'}
            frontfact=(r0.^2.*spl.*Nb.*sigh.*bx.*by)./(4.*sqrt(pi).*(relbeta.^2).*(relgamma.^4).*sigxb.^2.*sigyb.^2.*sigs.*sigp);
            
        otherwise
            frontfact=(r0.^2.*spl.*Nb)./(8.*pi.*(relgamma.^2).*sigs.*sqrt(...
                (sigx.^2).*(sigy.^2)-sigp.^4.*Dx.^2.*Dy.^2).*um).*2.*FpiWfact;
            
    end
    contributionsTL(:,dppcolnum)=frontfact.*val;
    
    L=zeros(length(lindata),1);
    for kk=1:length(lindata);L(kk)=lindata(kk).Length;end
    Tlcol(dppcolnum)=1/(1/sum(L)*sum(contributionsTL(:,dppcolnum).*L));
    
end

TLT=length(Tlcol)/(sum(1./Tlcol));