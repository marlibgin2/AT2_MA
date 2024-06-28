function ACHRO_SP=splitlat(ACHRO,split)
    splits=1/split*ones(1,split);
    ACHRO_SP={};
    for i=1:numel(ACHRO)
        if (ACHRO{i}.Length>0.0)
            element=atdivelem(ACHRO{i},splits);
            ACHRO_SP=[ACHRO_SP;element];
        else
            ACHRO_SP=[ACHRO_SP;ACHRO{i}];
        end
    end