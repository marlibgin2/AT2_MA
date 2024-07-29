function pp=ppool(c,n)
% Starts parpool in cluster c wtih n workers
% turns off warnings from all workers
% in aurora nworkers = n_nodes*24 (max n_npdes is 16)
% in cosmos nworkers = n_nodes*48 (up top 10 nodes tested so far)
pp=parpool(c,n);
pctRunOnAll warning off
end