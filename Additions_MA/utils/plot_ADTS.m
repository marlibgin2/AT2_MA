function diffu = plot_fma(x,y,qx,qy,dqx,dqy,daplotrange,tunerange,nturn,nxpoints,nypoints)

nx = length(unique(x));
ny = length(unique(y));

darange = [min(x) max(x) min(y) max(y)];

dnux=.5e-3*4;
dnuy=.5e-3*4;

xstep=(darange(2)-darange(1))/(nx-1);
ystep=(darange(4)-darange(3))/(ny-1);

% building the arrays to be plotted
for icounter = 1:length(qx)
        if (qx(icounter) ~= 0 || qy(icounter) ~= 0)
            % the small rectangule size may vary depending on dx dy
            dpx=xstep/2; dmx=xstep/2;
            dpy=ystep/2; dmy=ystep/2;
            
            % the rectangule is described in clockwise way from bottom left point
            x0pos(:,icounter) = ...
                [x(icounter)-dmx; x(icounter)-dmx; x(icounter)+dpx; x(icounter)+dpx];
            y0pos(:,icounter) = ...
                [y(icounter)-dmy; y(icounter)+dpy; y(icounter)+dpy; y(icounter)-dmy];
            
            % buiding rectangular area in the tune plot around each tune point
            nux=qx(icounter);
            nuy=qy(icounter);
            nuxpos(:,icounter) = ...
                [nux-dnux; nux-dnux; nux+dnux; nux+dnux];
            nuypos(:,icounter) = ...
                [nuy-dnuy; nuy+dnuy; nuy+dnuy; nuy-dnuy];
            
            
            % determining the diffusion coefficient
            delta_tune= dqx(icounter)^2+dqy(icounter)^2;
            if (delta_tune > 0)
                diffu(icounter) = sqrt(delta_tune);
                diffu(icounter) = log10(diffu(icounter)/(nturn/2));
%                 if (diffu(icounter) < (-10))
%                     diffu(icounter) = -10;
%                 end
            else
                diffu(icounter) = -10;
            end
            diffuvec(1:4,icounter) = diffu(icounter);
            diffuvec(1:4,icounter) = diffu(icounter);
        else
            x0pos(:,icounter)  = [100; 100; 100; 100];			
            y0pos(:,icounter)  = [100; 100; 100; 100];			
            nuxpos(:,icounter) = [0; 0; 0; 0];
            nuypos(:,icounter) = [0; 0; 0; 0];
            diffuvec(:,icounter) = [-10; -10; -10; -10];
            diffuvec(:,icounter) = [-10; -10; -10; -10];
        end
end

% list = find(diffuvec==-10);
% x0pos(list) = NaN;

% figure;
figure(1);
hold off
fill(x0pos*1e3, y0pos*1e3, diffuvec);
axis(daplotrange*1e3);
caxis([-10 -3]);
grid on;
hold on;
shading flat;
colormap('jet');
colorbar;
tstr1 = ['nturn = ', num2str(nturn)];
% title(tstr1);
xlabel('x (mm)');
ylabel('y (mm)');

figure(11); % nux vs x 
hold off
len = length(nuxpos);
nuxpos_r(1,1:len) = nuxpos(1,1:len);
nuxpos_r(2,1:len) = nuxpos(3,1:len);
nuxpos_r(3,1:len) = nuxpos(4,1:len);
nuxpos_r(4,1:len) = nuxpos(2,1:len);

fill(x0pos*1e3, nuxpos_r, diffuvec);
axis([-10 10 0. 0.5])

caxis([-10 -3]);
grid on;
hold on;
shading flat;
colormap('jet');
colorbar;
%tstr1 = ['nturn = ', num2str(nturn)];
% title(tstr1);
xlabel('x (mm)');
ylabel('\nu_x');

sely0 = abs(y)<0.1e-3; 
figure(1100)
plot(x(sely0),qx(sely0),'r.'); 
xlabel('x (mm)');
ylabel('\nu_x');
title('\delta\nu_x / \deltax')


figure(12); % nux vs y 
hold off
len = length(nuxpos);
nuxpos_r = nuxpos;
nuxpos_r = nuxpos;
nuxpos_r = nuxpos;
nuxpos_r = nuxpos;

fill(y0pos*1e3, nuxpos_r, diffuvec);
axis([0 3.5 0.4 0.5])

caxis([-10 -3]);
grid on;
hold on;
shading flat;
colormap('jet');
colorbar;
%tstr1 = ['nturn = ', num2str(nturn)];
% title(tstr1);
xlabel('y (mm)');
ylabel('\nu_x');

selx0 = abs(x)<0.1e-3; 
figure(1200)
plot(y(selx0),qx(selx0),'r.'); 
xlabel('y (mm)');
ylabel('\nu_x');
title('\delta\nu_x / \deltay')

% figure
figure(2)
hold off; 
fill(nuxpos, nuypos, diffuvec);
axis(tunerange);
caxis([-10 -3]);
hold on;
shading flat;
colormap('jet');
colorbar;
plot_net(5,tunerange(1),tunerange(2),tunerange(3),tunerange(4));

tstr2 = ['Frequency Map'];
% title(tstr2);
xlabel('\nu_x');
ylabel('\nu_y');

figure(22); % nux vs x 
hold off
len = length(nuypos);
nuypos_r(1,1:len) = nuypos(1,1:len);
nuypos_r(2,1:len) = nuypos(3,1:len);
nuypos_r(3,1:len) = nuypos(4,1:len);
nuypos_r(4,1:len) = nuypos(2,1:len);

fill(y0pos*1e3, nuypos_r, diffuvec);
axis([0 4.5 0.2 0.4])

caxis([-10 -3]);
grid on;
hold on;
shading flat;
colormap('jet');
colorbar;
%tstr1 = ['nturn = ', num2str(nturn)];
% title(tstr1);
xlabel('y (mm)');
ylabel('\nu_y');

selx0 = abs(x)<0.1e-3; 
figure(2200)
plot(y(selx0),qy(selx0),'r.'); 
xlabel('y (mm)');
ylabel('\nu_y');
title('\delta\nu_y / \deltay')

figure(21); % nuy vs x 
hold off
len = length(nuxpos);
nuypos_r = nuypos;
nuypos_r = nuypos;
nuypos_r = nuypos;
nuypos_r = nuypos;

fill(x0pos*1e3, nuypos_r, diffuvec);
axis([-10 10 0.2 0.4])

caxis([-10 -3]);
grid on;
hold on;
shading flat;
colormap('jet');
colorbar;
%tstr1 = ['nturn = ', num2str(nturn)];
% title(tstr1);
xlabel('x (mm)');
ylabel('\nu_y');

sely0 = abs(y)<0.1e-3; 
figure(2100)
plot(x(sely0),qy(sely0),'r.'); 
xlabel('x (mm)');
ylabel('\nu_y');
title('\delta\nu_y / \deltax')
