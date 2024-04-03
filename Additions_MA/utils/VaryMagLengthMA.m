function R=VaryMagLengthMA(R,DL,indP,indM)
% changes length of two mag families, keeping total length constant
ndp=length(indP);
ndm=length(indM);
R=setcellstruct(R,'Length',indP,getcellstruct(R,'Length',indP)+DL/ndp);
R=setcellstruct(R,'Length',indM,getcellstruct(R,'Length',indM)-DL/ndm);
return