function [C,~] = CONfcn(X, rrr)

twiss =  gettwiss(rrr, 0.0);
hi    = abs(twiss.eta(1));
C(1)  = hi-100e-6; 


end

