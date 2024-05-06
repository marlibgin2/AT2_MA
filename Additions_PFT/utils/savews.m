dt=sprintf('%s',datetime);
fn=strcat('WS_',dt);
save(fn,'-regexp','^(?!(c)$).');
