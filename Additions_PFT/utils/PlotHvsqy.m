figure; 
semilogy(RbTuneScan2.data(:,58),RbTuneScan2.data(:,J+23),'o');
xlabel('qy');
ylabel('h');
grid on; 
legend(hnames{J}); hold on;