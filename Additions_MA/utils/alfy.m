function [ay]=alfy(Seq,indx)
% get value of bety for  Seq(indx)

T=twissring(Seq,0,indx);
a=cat(1,T.alpha);
ay=a(:,2)';

end