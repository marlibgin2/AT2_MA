function plot_net_frac(order,qxmin,qxmax,qzmin,qzmax)
%PLOT_NET Plot the tune resonance diagram of desired order
% plot_net(order,qxmin,qxmax,qzmin,qzmax)
%
% This function plot the net of resonance up to order 'order' appearing in 
% the tune window (qxmin, qxmax, qzmin, qzmax). It plots the resonance net
% on top of an existing figure.
%
% INPUT:
%   order     maximum order of resonance net to plot
%   qxmin     minimum H tune (fractional)
%   qxmax     maximum H tune (fractional)
%   qzmin     minimum V tune (fractional)
%   qzmax     maximum V tune (fractional)
% OUTPUT:
%   the resonance net plot is superimposed to an existing plot
%
% Author: RB 01/2004
%

% colour code *** dull ***
cc='kkbrgmcyyyyyyyyyyyyyyyyyyyyyyyy';

% search for the straight lines crossing the input window
for ordercolor=order:-1:1
    for m=-ordercolor:ordercolor
        for n=-ordercolor:ordercolor
            for p=-ordercolor:ordercolor
                if(n == 0)
                    if(m ~= 0)
                        line([p/m,p/m], [0,1],'Color',cc(abs(m)));
                    end    
                else               
                    qy0=p/n;
                    qy1=(p-m)/n;
                    if(qy0*qy1 <= 0)
                        if((abs(m)+abs(n)) <= order)
                            line([0,1], [qy0,qy1],'Color',cc(abs(m)+abs(n)));
                        end    
                    else
                        if((qy0 >=0) & (qy0 <=1))
                            if((abs(m)+abs(n)) <= order)
                                line([0,1], [qy0,qy1],'Color',cc(abs(m)+abs(n)));
                            end
                        end
                        if((qy1 >=0) & (qy1 <=1))
                            if((abs(m)+abs(n)) <= order)
                                line([0,1], [qy0,qy1],'Color',cc(abs(m)+abs(n)));
                            end
                        end
                    end
                end
            end
        end 
    end
end
axis([qxmin,qxmax,qzmin,qzmax]);
