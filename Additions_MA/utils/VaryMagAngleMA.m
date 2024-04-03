%function R=VaryMagAngleMA(R,BendVal,indP,indM)
function R=VaryMagAngleMA(R,BendVal,famAB, famTGB, SumABTGB)
% changes angle of two mag families, keeping total bending constant
% typically used to alter Anti-Bends compensating with Trensverse Gradient
% bends
% ndp=length(indP);
% ndm=length(indM);
%ndp=1; ndm=1;
%R=setcellstruct(R,'BendingAngle',indP,getcellstruct(R,'BendingAngle',indP)+Dtheta/ndp);
%R=setcellstruct(R,'BendingAngle',indM,getcellstruct(R,'BendingAngle',indM)-Dtheta/ndm);

indAB  = findcells(R,'FamName',famAB);
indTGB = findcells(R,'FamName',famTGB); 

R=setcellstruct(R,'BendingAngle',indAB,BendVal*ones(size(indAB))); % set now AB to the new value
R=setcellstruct(R,'BendingAngle',indTGB,SumABTGB-BendVal); % remove DeltaAB from TGB


%DeltaBendVal = BendVal - getcellstruct(R,'BendingAngle',indP); % vary AB angle and compute DeltaAB
 

return