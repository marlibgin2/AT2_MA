function [lindata,avebeta,avemu,avedisp,nu,xsi]=atavedata(ring,dpp,refpts,varargin)
%ATAVEDATA       Average of optical functions on selected elements
%
%[LINDATA,AVEBETA,AVEMU,AVEDISP,TUNES,CHROMS]=ATAVEDATA(RING,DPP,REFPTS)
%
%LINDATA : Identical to ATLINOPT output
%AVEBEA :  Average Beta functions
%AVEMU :   Average phase advance
%AVEDISP : Average dispersion
%TUNES : Vector of tunes
%CHROMS : Vector of chromaticites
%
%[LINDATA,AVEBETA,AVEMU,AVEDISP,TUNES,CHROMS]=ATAVEDATA(RING,DPP,REFPTS,ORBITIN)
%    does not search for closed orbit. instead ORBITIN is used

lr=length(ring)+1;
if islogical(refpts)
    refs=[refpts(:);false(lr-length(refpts),1)];
else
    refs=false(lr,1); % lr
    refs(refpts)=true;
end
long=atgetcells(ring,'Length',@(elem,lg) lg>0) & refs(1:end-1); %lr-1
needed=refs | [false;long]; %lr
[lind,nu,xsi]=atlinopt(ring,dpp,find(needed),varargin{:}); %needed

lindata=lind(refs(needed)); %refpts
avebeta=cat(1,lindata.beta); %refpts
avemu=cat(1,lindata.mu); %refpts
avedisp=cat(2,lindata.Dispersion)'; %refpts

if any(long)
    initial=[long(needed(1:end-1));false]; %needed
    final=[false;initial(1:end-1)]; %needed

    lg=initial(refs(needed)); % refpts
    L=atgetfieldvalues(ring(long),'Length'); %long
    
    beta0=avebeta(lg,:); %long
    alpha0=cat(1,lind(initial).alpha); %long
    mu0=avemu(lg,:); %long
    disp0=avedisp(lg,:); %long
    
    beta1=cat(1,lind(final).beta); %long
    alpha1=cat(1,lind(final).alpha); %long
    mu1=cat(1,lind(final).mu); %long
    disp1=cat(2,lind(final).Dispersion)'; %long
    
    L2=[L L]; %long
    avebeta(lg,:)       = betadrift(beta0,beta1,alpha0,L2);
    avemu(lg,:)         = 0.5*(mu0+mu1);
    avedisp(lg,[1 3])   = (disp1(:,[1 3])+disp0(:,[1 3]))*0.5;
    avedisp(lg,[2 4])   = (disp1(:,[2 4])+disp0(:,[2 4]))*0.5;
    
    foc  = atgetcells(ring(long),'PolynomB',@(el,polb) length(polb)>=2 && polb(2)~=0); %long
    bnd  = atgetcells(ring(long),'BendingAngle',@(el,ang) ang ~= 0 );
    if any(foc)
        qp=false(size(long));
        qp(long)=foc;
        bend = false(size(long));
        bend(long) = bnd;
        K=zeros(size(L)); %long
        K(foc)=atgetfieldvalues(ring(qp),'PolynomB',{2});
        K2=[K -K]; %long
        Linvrho = zeros(size(L));
        Linvrho(bnd)=atgetfieldvalues(ring(bend),'BendingAngle');
        sel=false(size(avebeta,1),1); %refpts
        sel(lg)=foc;
        avebeta(sel,:)=betafoc(beta1(foc,:),alpha0(foc,:),alpha1(foc,:),K2(foc,:),L2(foc,:));
        avedisp(sel,:)=dispfoc(disp0(foc,:),disp1(foc,:),K2(foc,:),L2(foc,:),Linvrho(foc,:)./L2(foc,:));
    end
end

    function avebeta=betadrift(beta0,beta1,alpha0,L)
        gamma0=(alpha0.*alpha0+1)./beta0;
        avebeta=0.5*(beta0+beta1)-gamma0.*L.*L/6;
    end

    function avebeta=betafoc(beta1,alpha0,alpha1,K,L)
        gamma1=(alpha1.*alpha1+1)./beta1;
        avebeta=0.5*((gamma1+K.*beta1).*L+alpha1-alpha0)./K./L;
    end

    function avedisp=dispfoc(dispp0,dispp1,K,L,invrho)
        % Analytical solution based on quadrupole linear transfer matrix
        % Approaches drift solution in limit K -> 0
        % While <D> = (1/2)*( D'(0) - D'(L) )/(K*L) is analytically
        % identical and correct, it misbehaves for is more sensitive to small errors.
        avedisp = dispp0;

        k = (invrho.^2 + K);
%        k = K;
        avedisp(:,[1 3]) = 0.5*( (dispp0(:,[1 3]) + dispp1(:,[1 3])) .* sin(sqrt(k).*L)./(sqrt(k).*L) + (dispp0(:,[2 4]) - dispp1(:,[2 4]))./sqrt(k) .* (1 - cos(sqrt(k).*L))./(sqrt(k).*L) );
        avedisp(:,[2 4]) = 0.5*( (dispp0(:,[2 4]) + dispp1(:,[2 4])) .* sin(sqrt(k).*L)./(sqrt(k).*L) - (dispp0(:,[1 3]) - dispp1(:,[1 3]))./sqrt(k) .* (1 - cos(sqrt(k).*L))./(sqrt(k).*L) );

% avedisp(:,[1 3]) = ( (dispp0(:,[1 3])) .* sin(sqrt(k).*L)./(sqrt(k).*L) + (dispp0(:,[2 4]))./sqrt(k) .* (1 - cos(sqrt(k).*L))./(sqrt(k).*L) );
% avedisp(:,[1 3]) = ( (dispp1(:,[1 3])) .* sin(sqrt(k).*L)./(sqrt(k).*L) - (dispp1(:,[2 4]))./sqrt(k) .* (1 - cos(sqrt(k).*L))./(sqrt(k).*L) );

%         avedisp(:,[1 3]) = (1./L).*( dispp0(:,[1 3])./sqrt(K).*sin(sqrt(K).*L) + dispp0(:,[2 4])./K .* (1 - cos(sqrt(K).*L)));
%         avedisp(:,[2 4]) = (1./L).*( dispp0(:,[2 4])./sqrt(K).*sin(sqrt(K).*L) - dispp0(:,[1 3])./K .* (1 - cos(sqrt(K).*L)));

%         avedisp(:,[1 3]) = 0.5*(dispp0(:,[2 4]) - dispp1(:,[2 4]))./(k.*L);

    end

% Bloody hell, we'll need the combined gradient bend equations too...
% obviously that messes with the dispersion.
% More general approach perhaps? 


end
