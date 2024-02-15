function res = RBParamScan_AT2(dTheta0,dTheta1,Ntheta,dk0,dk1,Nk)
%RBParamScan Scan Rb angle and Focusing strength calculates beam parameters for
%   a ring composed only of unit cells

res.Ntheta = Ntheta;
res.Nk     = Nk;

for i=1:49
    NaNs(i)=NaN;
end
res.data=[];
formatSpec = 'i = %3d  j = %3d  dTheta = %8.3f mrad Krb = %12.7f m**-2 \n';
for i=0:Ntheta-1
    dTheta = dTheta0 +(dTheta1-dTheta0)/(Ntheta-1)*i;
    for j=0:Nk-1
        if (Nk>1)
            krb    = dk0 + (dk1-dk0)/(Nk-1)*j;
        else
            krb=dk0;
        end
                  
        pars=RBParam(dTheta,krb);
        %if (pars(12)>=3.0)
        %   pars(1:21)=NaN;
        %end
        res.data=cat(1,res.data,real(pars));
        fprintf(formatSpec,i,j,dTheta*1000,krb);
    end
end


