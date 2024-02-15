for i=1:121;
    habs(i,1)=RbTuneScan2.data(i,57);
    habs(i,2)=sum(RbTuneScan2.data(i,[27 28 31 32 33 39 40 41]));
end;
figure; semilogy(habs(:,1),habs(:,2),'o'),xlabel('qx');grid on;