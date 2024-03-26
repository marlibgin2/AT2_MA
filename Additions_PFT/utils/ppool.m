function ppool(c,n)
% Starts parpool in cluster c wtih 24*n workers
% turns off warnings from all workers
parpool(c,24*n);
pctRunOnAll warning off
end