function DFT_reg_z_interp(path, shiftpath, refchannel, scale, Nchunks, varargin)
    
    p = inputParser;
    addOptional(p,'edges',[0,0,0,0]);
    addOptional(p,'Nt',[]); %how many volumes to register (if empty, does all)
    addOptional(p,'optotune','true'); %whether to apply the optotune transformation
    addOptional(p,'tforms_optotune',[]);
    addOptional(p,'reftype','mean'); %options are 'median' or 'mean' for projecting the reference volume
    addOptional(p,'blurfactor',1); %width of gaussian to blur for DFT reg
    addOptional(p,'keepingfactor',0.95); %what proportion of frame to accountf or with shifts
    addOptional(p,'planescorr',3); %how many planes up and down to search for Z DFT reg
    addOptional(p,'save',true); %save the DFT and optotune shifts
    addOptional(p,'debug',false); %save intermediate tiffs
    
    parse(p,varargin{:});
    p = p.Results;

    info = pipe.io.sbxInfo(path);
    fdir = fileparts(path);
    
    Nx = info.sz(1) - p.edges(3) - p.edges(4);
    Ny = info.sz(2) - p.edges(1) - p.edges(2);
    Nz = size(info.otwave,2);
    
    if isempty(p.Nt)
        p.Nt = info.nframes / numel(info.otwave);
    end
        
    chunkframes = Nz*floor(p.Nt / Nchunks);

    RS = cell(3,Nchunks);
    CS = cell(3,Nchunks);
    ZS = cell(1,Nchunks);
    
    H = parfor_progressbar(Nchunks,'DFT registration');
    
    parfor chunk = 1:Nchunks
        %1) load reference chunk
        raw_chunk = pipe.imread(path,chunkframes*(chunk-1)+1, chunkframes,refchannel,[]);
        raw_chunk = raw_chunk(p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),:);
        raw_chunk = reshape(raw_chunk,Nx,Ny,Nz,[]);
        raw_chunk = imresize(raw_chunk,1/scale);
    
        %2) apply optotune registration, if desired
        if p.optotune
            unwarped_chunk = ApplyOptotuneWarp(raw_chunk,p.tforms_optotune);
        else
            unwarped_chunk = raw_chunk;
        end
        
        if p.debug
            pipe.io.writeTiff(unwarped_chunk,strcat(fdir,filesep,'unwarped_chunk.tif'));
            fprintf('written unwarped chunk\n');
        end
        
        %3) first round of XY DFT registration
        %make single target volume for entire chunk
        ref1 = defineReference(unwarped_chunk,size(unwarped_chunk,4),p.reftype);
        %calculate DFT shift between each volume and target vol
        [RS1, CS1] = DetermineXYShiftsFBS(unwarped_chunk,p.blurfactor,p.keepingfactor,ref1);
        %apply the shift (need it for subsequent steps)
        chunk_reg1 = ApplyXYShiftsFBS(unwarped_chunk,RS1,CS1);
        
        if p.debug
            pipe.io.writeTiff(chunk_reg1,strcat(fdir,filesep,'chunk_reg1.tif'));
            fprintf('written first XY DFT\n');
        end
        
        %clear raw chunk to save space
        unwarped_chunk = [];
        
        %4) second round of Z interpolation registration
        %make new reference volume based on XY DFT registration (3)
        ref2 = defineReference(chunk_reg1,size(chunk_reg1,4),p.reftype);
        %calculate interpolated zshifts by fitting polynomial 
        [RS2,CS2,ZS1] = ComputeZshiftInterpolateFBS(ref2, chunk_reg1, p.planescorr, p.edges);
        %apply zshift
        chunk_reg2 = ApplyZShiftInterpolateFBS(chunk_reg1, ZS1, CS2, RS2);
        if p.debug
            pipe.io.writeTiff(chunk_reg2,strcat(fdir,filesep,'chunk_reg2.tif'));
            fprintf('written XYZ DFT\n');
        end
        %find NaNs and replace with zeros
        nan_idx = find(isnan(chunk_reg2));
        chunk_reg2(nan_idx) = 0;
        %clear registration from step 1
        chunk_reg1 = [];
        
        %5) third round of XY DFT registration
        %make new reference volume based on XYZ DFT (4)
        ref3 = defineReference(chunk_reg2,size(chunk_reg2,4),p.reftype);
        %calculate DFT shift between each volume and target vol
        [RS3, CS3] = DetermineXYShiftsFBS(chunk_reg2,p.blurfactor,p.keepingfactor,ref3);
           
        if p.debug
            chunk_reg3 = ApplyXYShiftsFBS(chunk_reg2,RS3,CS3);
            pipe.io.writeTiff(chunk_reg3,strcat(fdir,filesep,'chunk_reg3.tif'));
            fprintf('written second XY DFT\n');
        end
        %clear raw chunk to save space
        
        
        chunk_reg2 = [];
        
        %6) combine Row and Column Shifts from (3), (4), and (5)
        RS(:,chunk) = {RS1,RS2,RS3};
        CS(:,chunk) = {CS1,CS2,CS3};
        ZS(chunk) = {ZS1};
        
        %7) save reference files for stitching later
        ref_all{chunk} = ref3;
        
        %iterate 
        H.iterate(1);
    end
    %Fix intra-chunk discontinuities
    ref_final = ref_all{1};
    %concatenate all of the reference volumes from above
    ref_cat = ref_all{1};
    for i = 2:Nchunks
        ref_cat = cat(4,ref_cat,ref_all{i});
    end
    
    %Compute XYZ shifts from each chunk's reference volume
    [RS_chunk,CS_chunk,ZS_chunk] = ComputeZshiftInterpolateFBS(ref_final, ref_cat, p.planescorr, p.edges);
    ref_cat = [];
    
    %convert local shift correction cells to matrix form
    CS1 = [CS{1,:}];
    CS2 = [CS{2,:}];
    CS3 = [CS{3,:}];
    RS1 = [RS{1,:}];
    RS2 = [RS{2,:}];
    RS3 = [RS{3,:}];
    ZS1 = [ZS{:}];

    %stretch the intra-chunk corrections to apply to every frame
    RS_chunk = imresize(RS_chunk,size(RS1),'nearest');
    CS_chunk = imresize(CS_chunk,size(CS1),'nearest');
    ZS_chunk = imresize(ZS_chunk,size(ZS1),'nearest');
    
    if p.debug
        tic
        raw_chunk = pipe.imread(path,chunkframes*(1-1)+1, chunkframes,refchannel,[]);
        ch = ApplyXYShiftsFBS(raw_chunk,RS1(:)*scale,CS1(:)*scale);
        fprintf('one shift applied\n');
        ch = ApplyZShiftInterpolateFBS(ch,ZS1(:),CS2(:)*scale,RS2(:)*scale);
        fprintf('two shifts applied\n');
        ch = ApplyXYShiftsFBS(ch,RS3(:)*scale,CS3(:)*scale);
        fprintf('three shifts applied\n');
        ch = ApplyXYShiftsFBS(ch,RS_chunk(:)*scale,CS_chunk(:)*scale);
        fprintf('four shifts applied\n');
        toc
        pipe.io.writeTiff(ch,strcat(fdir,filesep,'chunk_reg_piecewise.tif'));
        fprintf('written final DFT\n');
    end
    
    %scale the optotune transforms
    if p.optotune
        tforms_optotune_full = p.tforms_optotune;
        for i = 1:size(p.tforms_optotune,2)
            tforms_optotune_full(i).T(3,1:2) = tforms_optotune_full(i).T(3,1:2)*scale;
        end
    else
        tforms_optotune_full = [];
    end
    
    %save the DFT and optotune transformations
    if p.save
%         save(strcat(fdir,filesep,'DFTShifts.mat'),'RS_final','CS_final','ZS_final','tforms_optotune_full');
        save(shiftpath,'RS1','CS1','ZS1','RS2','CS2','RS3','CS3',...
            'RS_chunk','CS_chunk','ZS_chunk','scale','tforms_optotune_full','-mat');
    end
    close(H);
end























