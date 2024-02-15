for i=1:121;
    habs(i,1)=RbTuneScan2.data(i,57);
    habs(i,2)=sum(RbTuneScan2.data(i,26:46));
end;
figure; plot(habs(:,1),habs(:,2),'-o'),xlabel('qy');