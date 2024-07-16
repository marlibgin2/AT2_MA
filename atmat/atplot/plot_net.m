function plot_net(order,qxmin,qxmax,qzmin,qzmax)
% It plot the tune resonance diagram upto desired order
% plot_net_new(order,qxmin,qxmax,qzmin,qzmax)
%% INPUT:
%   order     maximum order of resonance net to plot
%   qxmin     minimum H tune
%   qxmax     maximum H tune
%   qzmin     minimum V tune
%   qzmax     maximum V tune
%% Usage examples
% plot_net(3,55.0,55.5,16.0,16.5);

%% History
% SKJ 2024/07/15 

nux =[qxmin,qxmax];
nuy = [qzmin,qzmax];
P=1;   % Super period
t=ceil(order*max([nux,nuy])/P);

%color code for resonance lines
% 1st order: cyan
% 2nd order : green
% 3rd order: red
% 4th order : blue
% 5th order: black
% 6th order : magenta

mm='cgrbkm';

hold on;
for o=order:-1:1
    ord=-o:o;
    k=1;
    lm=[];
    for i=1:length(ord)
        l=ord(i);
        for j=1:length(ord)
            m=ord(j);
            if (abs(l)+abs(m))==o
                lm(k,:)=[l,m];
                k=k+1;
            end
        end
    end
    for ii=1:max(size(lm))
        if lm(ii,2)~=0
            for z=0:t
                for jj=1:length(nux)
                    nuy1(jj)=(z*P-lm(ii,1).*nux(jj))/lm(ii,2);
                end
                plot(nux,nuy1,mm(o), 'LineWidth',1)
            end
        end
        if lm(ii,2)==0
            for z=0:t
                for kk=1:length(nuy)
                    nux1(kk)=z*P/lm(ii,1);
                end
                plot(nux1,nuy,mm(o),'LineWidth',1)
            end
        end
    end
end
 ylim(nuy);
 xlim(nux); 

