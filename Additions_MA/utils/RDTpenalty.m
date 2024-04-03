function rdtp = RDTpenalty(RDT)

nfield = length(fieldnames(RDT));

rdtp = 0; fn = fieldnames(RDT);
for i = 1:nfield
    ff = fn{i};
    
    if  ~strcmpi(ff,'h11001') && ~strcmpi(ff,'h00111') % not the chromaticity! 
        rdtp = rdtp + abs(RDT(end).(ff));
    end

end

