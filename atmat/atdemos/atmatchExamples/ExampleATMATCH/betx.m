function [bx]=betx(Seq,indx)
% get value of betx for  Seq(indx)

bx=1e+19;
try
    T=twissring(Seq,0,indx);
    b=cat(1,T.beta);
    bx=b(:,1)';
catch
    disp('betax twiss results is unphysical ...')
end
end