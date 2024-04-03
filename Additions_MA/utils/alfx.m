function [ax]=alfx(Seq,indx)
% get value of betx for  Seq(indx)

T=twissring(Seq,0,indx);
a=cat(1,T.alpha);
ax=a(:,1)';

end