function qxvsqy = ucXvsYTune(LatticeOptData,qy0,dqy,Nqy,...
                             emittance,nit,TolTunes,verbosef)
% Generates table of horizontal unit cell tune need to achieve a certain
% emittance as a fucntion of vertical unit cell tune
%   

k=1;
qx=NaN(2*Nqy+1,1);
qy=NaN(2*Nqy+1,1);
for j=-Nqy:Nqy
        qy(k) = qy0 + j*dqy;
        fprintf ('qy = %9.6f \n',qy(k) );
        qx(k)= uctunexvsemit(LatticeOptData,qy(k),emittance,...
                      nit,TolTunes,verbosef);
        k=k+1;
end
qxvsqy.qy0 = qy0;
qxvsqy.dqy = dqy;
qxvsqy.emittance = emittance;
qxvsqy.nit = nit;
qxvsqy.TolTunes = TolTunes;
qxvsqy.LatticeOptData=LatticeOptData;
qxvsqy.data = cat(2,qy,qx);
end

