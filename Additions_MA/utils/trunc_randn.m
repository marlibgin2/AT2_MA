function Xout = trunc_randn(Nel0, Nsigma)
%
% generate a normal gaussian distribution truncated at Nsigma
%
if nargin<2
    Nsigma=2;
end
Nel = Nel0;
Xout = [];
while Nel>0
    X = randn(Nel,1);
    X = X(abs(X)<Nsigma);
    Xout = [Xout; X]; %#ok<AGROW> 
    Nel   = Nel - max(size(X));
end

end

