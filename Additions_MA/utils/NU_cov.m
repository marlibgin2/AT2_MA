function [covM, h, Ae] = NU_cov(nturn) 

clear Mfx Mfy Mdfx Mdfy fx fy dfx dfy 
fmap = 'fma.out.out';
[x0 y0 fx fy dfx dfy] = textread(fmap,'%f %f %f %f %f %f','headerlines',1);
%Mfx  = reshape(fx,160,160);
%Mfy  = reshape(fy,160,160);
%Mfx  = Mfx(dfx~=0);
%Mfy  = Mfy(dfx~=0);
dnu = sqrt(dfx.^2+dfy.^2);
dnu = log10(dnu/1056);
%Mdnu = reshape(-dnu,160,160);
%Mdfx = reshape(dfx,160,160);
%Mdfy = reshape(dfy,160,160);
%Mdfx = Mdfx(dfx~=0);
%Mdfy = Mdfy(dfx~=0);

%%%%fx=fx(dfx~=0); fx=fx(~isnan(fx)); fx=fx(~isnan(fy))
%%%%fy=fy(dfx~=0); fy=fy(~isnan(fx)); fy=fy(~isnan(fy)) 

all=(1:length(fx))';
ax = find(fx==0|isnan(fx));
ay = find(fy==0|isnan(fy));
bx = find(dfx==0|isnan(dfx));
by = find(dfy==0|isnan(dfy));

a = union(ax,ay); 
b = union(bx,by);
ab = union(a,b);
g=setdiff(all,ab);

fx=fx(g); fy=fy(g); dfx=dfx(g); dfy=dfy(g); 

%dfx=dfx(dfx~=0); dfx=dfx(~isnan(fx));
%dfy=dfy(dfx~=0); dfy=dfy(~isnan(fx));
dnu = sqrt(dfx.^2+dfy.^2);
dnu = log10(dnu/nturn);
dnu=-dnu;
dnu=dnu/sum(dnu); % fractionary weight
w_nux  = fx.*dnu;
w_nuy  = fy.*dnu;
w_nux2 = (fx.^2).*dnu;
w_nuy2 = (fy.^2).*dnu;
w_nuxy = (fx.*fy).*dnu;

% nux    = fx;
% nuy    = fy; 
% nux2   = fx.*fx;
% nuy2   = fy.*fy;

Sxx    = sum(w_nux2) - (sum(w_nux))^2;
Syy    = sum(w_nuy2) - (sum(w_nuy))^2;
Sxy    = sum(w_nuxy) - sum(w_nux)*sum(w_nuy);

covM = [Sxx Sxy; ...
       Sxy Syy];
Ae = pi*sqrt(Sxx * Syy);  
h  = error_ellipse(covM,[sum(w_nux), sum(w_nuy)]);
end



