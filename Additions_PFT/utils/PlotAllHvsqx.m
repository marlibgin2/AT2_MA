figure; 
semilogy(RbTuneScan2.data(:,57),RbTuneScan2.data(:,26),'o');
xlabel('qx');ylabel('h');grid on; 
legend(hnames{3}); hold on;

for i=28:46
    semilogy(RbTuneScan2.data(:,57),RbTuneScan2.data(:,i),'o');
    legend(hnames{i-28+4});
end;

