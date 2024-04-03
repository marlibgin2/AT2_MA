function NEWRING = ataddshift(RING,index_list,dx,dy)

for i = 1:numel(index_list)
    dx0(i) = -RING{index_list(i)}.T1(1);
    dy0(i) = -RING{index_list(i)}.T1(3);
end
S  = size(dx);
S0 = size(dx0);
if S(1)==S0(2)
    dx0=dx0';
    dy0=dy0';
end

NEWRING = atsetshift(RING, index_list, dx+dx0, dy+dy0);
