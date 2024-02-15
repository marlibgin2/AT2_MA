Xmax = 0.015;
Ymax = 0.005;
npdax = 32;
npday = 32;

dx = Xmax/npdax;
dy = Ymax/npday;

np = (2*npdax+1)*(npday+1);
X0da = zeros(np,1);
Y0da = zeros(np,1);

k= 1;
for i=0:npday
    for j=0:2*npdax
        X0da(k) = -Xmax+dx*j;
        Y0da(k) = dy*i;
        k=k+1;
    end
end



