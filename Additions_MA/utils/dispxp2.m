function [dxp]=dispxp2(Seq,indx)
% get value of horizontal dispersion for  Seq(indx)

% try to compute dispersion, if it fails give it an absurd value
try
    [~,dxp,~,~]=getDispersion2(Seq,indx);
catch
    dxp=ones(1,length(indx))*1e19;
end

end