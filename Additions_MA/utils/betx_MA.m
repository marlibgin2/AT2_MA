function [bx]=betx_MA(Seq,indx)
% get value of betx for  Seq(indx)

% try to compute beta, if it fails give it an absurd value
try
    T=twissring(Seq,0,indx);
    b=cat(1,T.beta);
    bx=b(:,1)';
catch
    disp('exception in betx found ... ')
    bx=ones(1,length(indx))*1e19;
end

end