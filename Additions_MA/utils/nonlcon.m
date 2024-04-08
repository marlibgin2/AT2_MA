function valore = nonlcon(X0)

h=ellipse(X0(1),X0(2),0,X0(3),0,'b',100,0);
xe=h.XData(h.YData>=0);
ye=h.YData(h.YData>=0);

DAtest  =  [0.0051    0; ...
    0.0048    0.0008; ...
    0.0044    0.0014; ...
    0.0050    0.0025; ...
    0.0043    0.0031; ...
    0.0034    0.0034; ...
    0.0028    0.0039; ...
    0.0024    0.0046; ...
    0.0016    0.0051; ...
    0.0008    0.0051; ...
    0.0000    0.0053; ...
   -0.0008    0.0052; ...
   -0.0017    0.0052; ...
   -0.0026    0.0051; ...
   -0.0030    0.0041; ...
   -0.0034    0.0034; ...
   -0.0037    0.0027; ...
   -0.0047    0.0024; ...
   -0.0071    0.0023; ...
   -0.0085    0.0013; ...
   -0.0107    0.0000];

DAtest = 1e3*DAtest;
DAtest(DAtest(:,2)>2)=2;

[in,on] = inpolygon(xe,ye,DAtest(:,1),DAtest(:,2));

if (~isempty(find(in==0)))
    valore = 1;
else
    valore=-1;
end



