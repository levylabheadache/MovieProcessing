function regref = DFT_warp_2D(path, shiftpath, refchannel, chunksize, varargin)
    p = inputParser;
    addOptional(p,'edges',[0,0,0,0]);
    addOptional(p,'Nt',[]); %how many volumes to register (if empty, does all)
    addOptional(p,'reftype','mean'); %options are 'median' or 'mean' for projecting the reference volume
    addOptional(p,'blurfactor',1); %width of gaussian to blur for DFT reg
    addOptional(p,'keepingfactor',0.95); %what proportion of frame to accountf or with shifts
    addOptional(p,'planescorr',3); %how many planes up and down to search for Z DFT reg
    addOptional(p,'save',true); %save the DFT and optotune shifts
    addOptional(p,'iterations',3);
    
    parse(p,varargin{:});
    p = p.Results;

    info = pipe.io.sbxInfo(path);
    fdir = fileparts(path);
    
    Nchan = info.nchan;
    Nx = info.sz(1) - p.edges(3) - p.edges(4);
    Ny = info.sz(2) - p.edges(1) - p.edges(2);
    
    if isempty(p.Nt)
        p.Nt = info.nframes / numel(info.otwave);
    end
    
    Nchunks = ceil(p.Nt / chunksize);
%     chunksize = floor(p.Nt / Nchunks);
    
    chunkidx = 1:chunksize:p.Nt;
    
    RS = cell(1,Nchunks);
    CS = cell(1,Nchunks);
    ref = cell(1,Nchunks);
    
    H = parfor_progressbar(Nchunks,'DFT registration');
    parfor j = 1:Nchunks
        
        rawchunk = zeros(Nchan,Nx,Ny,chunksize);
        rawchunk = pipe.io.sbxRead(path,chunkidx(j),chunksize, -1, []);

        %initialize ref image
        if Nchan == 2
            regtemp = rescale(squeeze(rawchunk(refchannel,:,:,:)));
        elseif Nchan == 1
            regtemp = rescale(rawchunk);
        end
        
        
        chunkref = mean(regtemp,3);

        Rtot = zeros(size(rawchunk,ndims(rawchunk)),1);
        Ctot = zeros(size(rawchunk,ndims(rawchunk)),1);

        %do DFT reg the number of iterations defined above (default 3)
        for i = 1:p.iterations
            refchunk = mean(regtemp,3);
            
            if i ~= p.iterations
                [Rtemp,Ctemp,regtemp] = DFT_reg(regtemp,chunkref,10);
            else 
                [Rtemp,Ctemp,~] = DFT_reg(regtemp,chunkref,10);
            end

            Rtot = Rtot+Rtemp;
            Ctot = Ctot+Ctemp;
            chunkref = mean(regtemp,3);
        end
        RS{j} = Rtot;
        CS{j} = Ctot;
        ref{j} = chunkref;
        H.iterate(1);
    end
    close(H);
    refcat = cat(ndims(ref{1})+1,ref{:});
    [Rref,Cref,regref] = DFT_reg(refcat,refcat(:,:,1),10);
    if p.save
        save(shiftpath,'RS','CS','Rref','Cref','ref');
    end
end

    
    
    
    