%dt=sprintf('%s',datetime);
%fn=strcat('/home/pedtav/Documents/Max4U/WorkspaceBU/WS_',dt);
fn=strcat('/home/pedtav/Documents/Max4U/WorkspaceBU/WS_',datestr(now,30));
save(fn,'-regexp','^(?!(c)$).');
