function [by]=bety(Seq,indx)
% get value of bety for  Seq(indx)

by=1e+19;
try
    T=twissring(Seq,0,indx);
b=cat(1,T.beta);
by=b(:,2)';
catch
    disp('betay twiss results is unphysical ...')
end
end