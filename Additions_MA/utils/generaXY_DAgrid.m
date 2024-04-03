function [DA, mDA] = generaXY_DAgrid(Rin, nA)
if nargin<2
    nA = 1;
end
Nx   = 50*2;
xmin = -12e-3; xmax = 12e-3; 
Ny   = 32*2; 
ymin =   0e-3; ymax = 7e-3; 

X0 = linspace(xmin, xmax, Nx)';
Y0 = linspace(ymin, ymax, Ny)';

[XX0, YY0] = meshgrid(X0, Y0);

%XX0 = repmat(X0,Ny,1);
%YY0 = vertcat(repmat(Y0,Nx,1),repmat(1e-3,25,1)

y0 = reshape(YY0',Nx*Ny,1);
x0 = reshape(XX0',Nx*Ny,1);

%[nup, xip] = check_main_params(Rin,nA); 
RP = check_main_params(Rin,nA); 
nup=RP.nup; xip=RP.xip;

nturns = 2000;
DA = calcDA_grid(Rin, x0, y0, nturns*nA, 0);

if ~isempty(DA)
mDA=(reshape(DA,Nx, Ny));
figure(336); clf; hold on 

h = pcolor(XX0',YY0',mDA); set(h, 'EdgeColor', 'none'); box on
axis([-12e-3 12e-3 0 7e-3])
xlabel('X (m)'); ylabel('Y (m)'); title(['DA - long straight (' num2str(nturns) 't)'])
text(-11.5e-3, 6e-3,['nu = (' num2str(nup(1),5) ',' num2str(nup(2),5) ')'],'color','w')
text(-11.5e-3, 5.5e-3,['xi = (' num2str(xip(1),5) ',' num2str(xip(2),5) ')'],'color','w')

end
