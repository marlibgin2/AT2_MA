%function R=VaryMagAngleMA(R,BendVal,indP,indM)
function R=VaryDIPGrad3(R,kDIP6,famDIP)
% changes angle of two mag families, keeping total bending constant
% typically used to alter Anti-Bends compensating with Trensverse Gradient
% bends
% ndp=length(indP);
% ndm=length(indM);
%ndp=1; ndm=1;
%R=setcellstruct(R,'BendingAngle',indP,getcellstruct(R,'BendingAngle',indP)+Dtheta/ndp);
%R=setcellstruct(R,'BendingAngle',indM,getcellstruct(R,'BendingAngle',indM)-Dtheta/ndm);
frac = [-0.000128 -0.000128 0.011759 -0.551829 -0.866059 -0.864858]/-0.864858;
frac = [frac, flip(frac)]'; frac=repmat(frac,20,1);
indDIP  = findcells(R,'FamName',famDIP); 
for i = 1:length(indDIP)
    K(i)       =  kDIP6 * frac(i); 
    %R{indDIP(i)}.K *GF;
end

R=setcellstruct(R,'PolynomB',indDIP,K,1,2); % set now AB to the new value
R=setcellstruct(R,'K',indDIP,K); % set now AB to the new value


return