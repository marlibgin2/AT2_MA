function R=CompAngMA(R,Dtheta,indP,indM)
% this is a CONSTRAINT function: when the 1st angel is increased the second
% decreases
% changes angle of two mag families, keeping total bending constant
% typically used to alter Anti-Bends compensating with Trensverse Gradient
% bends
ndp=length(indP);
ndm=length(indM);
R=setcellstruct(R,'BendingAngle',indP,getcellstruct(R,'BendingAngle',indP)+Dtheta/ndp);
R=setcellstruct(R,'BendingAngle',indM,getcellstruct(R,'BendingAngle',indM)-Dtheta/ndm);
return