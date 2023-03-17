function tifname = MakeTifName(base,Nchan,Nz,Nt,c,z,t)
    tifname = base;
    
    if Nchan > 1 || Nz > 1 || Nt > 1
        tifname = strcat(tifname,'_');
    
    if Nchan > 1
        tifname = strcat(tifname,'C',sprintf('%03d',c));
    end
    
    if Nz > 1
        tifname = strcat(tifname,'Z',sprintf('%03d',z));
    end
    
    if Nt > 1
        tifname = strcat(tifname,'T',sprintf('%03d',t));
    end
    
end