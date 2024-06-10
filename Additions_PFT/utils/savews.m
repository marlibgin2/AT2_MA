dt=sprintf('%s',datetime);
fn=strcat('/home/pedtav/Documents/Max4U/WorkspaceBU/WS_',dt);
save(fn,'-regexp','^(?!(c)$).');
