function [by]=bety_MA(Seq,indx)
% get value of bety for  Seq(indx)

% try to compute beta, if it fails give it an absurd value
try
    T=twissring(Seq,0,indx);
    b=cat(1,T.beta);
    by=b(:,2)';
catch
    disp('exception in bety found ... ')
    by=ones(1,length(indx))*1e19;
end

end