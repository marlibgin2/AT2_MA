%function R=VaryMagAngleMA(R,BendVal,indP,indM)
function R=VaryDIPGrad(R,GradFac,kDIPstart,famDIP)
% changes angle of two mag families, keeping total bending constant
% typically used to alter Anti-Bends compensating with Trensverse Gradient
% bends
% ndp=length(indP);
% ndm=length(indM);
%ndp=1; ndm=1;
%R=setcellstruct(R,'BendingAngle',indP,getcellstruct(R,'BendingAngle',indP)+Dtheta/ndp);
%R=setcellstruct(R,'BendingAngle',indM,getcellstruct(R,'BendingAngle',indM)-Dtheta/ndm);

indDIP  = findcells(R,'FamName',famDIP); 

R=setcellstruct(R,'PolynomB',indDIP,GradFac*kDIPstart,1,2); % set now AB to the new value


%DeltaBendVal = BendVal - getcellstruct(R,'BendingAngle',indP); % vary AB angle and compute DeltaAB
 

return