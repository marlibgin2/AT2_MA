function [dx]=dispx2(Seq,indx)
% get value of horizontal dispersion for  Seq(indx)

% try to compute dispersion, if it fails give it an absurd value
try
    [dx,~,~,~]=getDispersion2(Seq,indx);
catch
    disp('exception in dispx2 found ... ')
    dx=ones(1,length(indx))*1e19;
end

end