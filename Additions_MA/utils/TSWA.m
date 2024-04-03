function TSWA(ring)

%
% MA - attempt at calculating ADTS follwing the recipe in 
% K. Soutome, H. Tanaka   PRAB20 064001 (2017)
%

nux_tilde = nux + 2*cxx.*Jx +   cxy.*Jy + 3*cxxx.*Jx.^2 + 2*cxxy.*Jx.*Jy;
nuy_tilde = nuy +   cxy.*Jx + 2*cyy.*Jy + cxxy.*Jx.^2;


end

function cxx
end