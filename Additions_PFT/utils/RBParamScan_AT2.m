function res = RBParamScan_AT2(dTheta0,dTheta1,Ntheta,dk0,dk1,Nk,chroms,nuc)
%RBParamScan Scan Rb angle and Focusing strength and calculates beam parameters for
%   a ring composed only of unit cells. Chromaticity is fit to the desired
%   vaue per cell. chroms is an array with two values for the horizontal
%   and vertical chromaticty per cell
%   here we 

res.Ntheta = Ntheta;
res.Nk     = Nk;
res.nuc    = nuc;

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
                  
        pars=RBParam_AT2(dTheta,nuc,chroms,krb);
        if (pars(12)>=3.0) % Checks for longitudinally unstable motion (Jx>3)
            pars(3:52)=NaN;
        end
        res.data=cat(1,res.data,real(pars));
        fprintf(formatSpec,i,j,dTheta*1000,krb);
    end
end


