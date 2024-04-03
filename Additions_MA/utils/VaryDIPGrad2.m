%function R=VaryMagAngleMA(R,BendVal,indP,indM)
function R=VaryDIPGrad2(R,GF,famDIP)
% changes angle of two mag families, keeping total bending constant
% typically used to alter Anti-Bends compensating with Trensverse Gradient
% bends
% ndp=length(indP);
% ndm=length(indM);
%ndp=1; ndm=1;
%R=setcellstruct(R,'BendingAngle',indP,getcellstruct(R,'BendingAngle',indP)+Dtheta/ndp);
%R=setcellstruct(R,'BendingAngle',indM,getcellstruct(R,'BendingAngle',indM)-Dtheta/ndm);

indDIP  = findcells(R,'FamName',famDIP); 
r       = kDIP6/kDIPstart0(6);

% % % for i=1:12
% % %     iDIP{i} = indDIP(linspace(i,1200-12+i,1200/12));
% % % end
% % % 
% % % for i = 1:12
% % %     R=setcellstruct(R,'PolynomB',iDIP{i},kDIPstart0(i).*r,1,2); % set now AB to the new value
% % % end

    
R=setcellstruct(R,'PolynomB',indDIP,kDIPstart0*r,1,2); % set now AB to the new value


return