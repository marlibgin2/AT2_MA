function elem=atepukick(fname,varargin)

% MIKm = epukick('MIK', nx, ny, XMIK, YMIK, PXMIK, PYMIK, 'ThinEPUPass');

[rsrc,nx,ny,XEPU,YEPU,PXEPU,PYEPU,method]      = decodeatargsMA({0,0,0,0,0,0,'ThinEPUPass'},varargin);

[nx,rsrc]          = getoption(rsrc,'nx',nx);
[ny,rsrc]          = getoption(rsrc,'nx',ny);
[XEPU,rsrc]        = getoption(rsrc,'XEPU',XEPU);
[YEPU,rsrc]        = getoption(rsrc,'YEPU',YEPU);
[PXEPU,rsrc]       = getoption(rsrc,'PXMIK',PXEPU);
[PYEPU,rsrc]       = getoption(rsrc,'PYMIK',PYEPU);

% Build the element
elem=atbaselem(fname,method,'Class','InjKick','Length',0, ...
    'NumX',nx,'NumY',ny,...
    'XGrid',  XEPU,  'YGrid',  YEPU, ...
    'PxGrid', PXEPU, 'PyGrid', PYEPU, ...
    'R1',diag(ones(6,1)),'R2',diag(ones(6,1)),...
    'T1',zeros(1,6),'T2',zeros(1,6),...
    rsrc{:});
end
