function CorrectExpt(expt, sbxInfo, regParam, varargin) % , tforms_optotune , Nchunks shiftPath,
IP = inputParser;
addRequired(IP, 'expt', @isstruct )
addRequired(IP, 'sbxInfo', @isstruct )
addRequired(IP, 'regParam', @isstruct )
%addParameter(IP, 'sbx', 'opt', @ischar) % which sbx file to use as the input
addParameter(IP, 'chunkSize', 15, @isnumeric)
addParameter(IP, 'blur', 1, @isnumeric); %width of gaussian to blur for DFT reg
addParameter(IP, 'keep', 0.95, @isnumeric); %what proportion of frame to account or with shifts
addParameter(IP, 'anchor',0, @isnumeric); %
addParameter(IP, 'zprojPath', '', @ischar)
%addParameter(IP, 'proj', true, @islogical);
addParameter(IP, 'range', [0.25, 0.85], @isnumeric);
parse(IP, expt, sbxInfo, regParam, varargin{:}); % tforms_optotune,  shiftPath,
%sbx_type = IP.Results.sbx;
chunkSize = IP.Results.chunkSize;
anchor = IP.Results.anchor;
blurFactor = IP.Results.blur;
keepFactor = IP.Results.keep;
%projToggle = IP.Results.proj;
proj_range = IP.Results.range;

if expt.Nplane > 1
    sbxInputPath = expt.sbx.opt;
else
    sbxInputPath = expt.sbx.cat;
end
if ~exist(sbxInputPath, 'file'), error('%s does not exist!', sbxInputPath); end
[refPMT, ~] = DeterminePMT(regParam.refChan, sbxInfo);
[usePMT, ~] = DeterminePMT('both', sbxInfo);
Npmt = numel(usePMT);
Nx = sbxInfo.sz(1); Ny = sbxInfo.sz(2);
Nx_crop = Nx-regParam.edges(3)-regParam.edges(4);
Ny_crop = Ny-regParam.edges(1)-regParam.edges(2);
resize_factor = 1/regParam.binXY;
rectify_start = round(sbxInfo.Nplane/2);
zProjRange = [1, 1];
if sbxInfo.Nplane > 1
    zProjRange = round(proj_range(1)*sbxInfo.Nplane):round(proj_range(2)*sbxInfo.Nplane);
end
zprojPath = sprintf('%s%s_zproj_%i-%i.tif', expt.dir, expt.name, zProjRange(1), zProjRange(end));
shiftPath = strcat(expt.dir, expt.name, '_dftshifts.mat');
if ~exist(shiftPath, 'file')
    [chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.totScan, 'size',chunkSize, 'allowPartial',true );
    RS = cell(1,Nchunk);  CS = cell(1,Nchunk);  ZS = cell(1,Nchunk);
    tic;
    w = waitbar(0,'DFT registration'); %H = parfor_progressbar(Nchunk,'DFT registration'); %
    for c = 1:Nchunk
        % load chunk of optotune-corrected data
        raw_chunk = readSBX(sbxInputPath, sbxInfo, chunkLims(c,1), chunkLength(c), refPMT ); % readSBX(path, info, k, N, pmt, optolevel)
        raw_chunk = raw_chunk(regParam.edges(3)+1:end-regParam.edges(4), regParam.edges(1)+1:end-regParam.edges(2), :);
        raw_chunk = reshape(raw_chunk, Nx_crop, Ny_crop, sbxInfo.Nplane, []);
        raw_chunk = imresize(raw_chunk, resize_factor);

        % rectify each volume with dft
        RS0 = zeros(sbxInfo.Nplane,chunkLength(c)); CS0 = zeros(sbxInfo.Nplane,chunkLength(c));
        if sbxInfo.Nplane > 1
            chunk_reg0 = nan(size(raw_chunk));
            parfor i = 1:chunkLength(c)
                [RS0(:,i), CS0(:,i), chunk_reg0(:,:,:,i)] = RectifyStack( raw_chunk(:,:,:,i), rectify_start) ; % DFT_rect( raw_chunk(:,:,:,i), , 4);
            end
        else
            chunk_reg0 = raw_chunk;
        end

        % first round of XY DFT registration
        ref1 = mean(chunk_reg0, 4); %defineReference(chunk_reg0, chunkSize, refType); % make single reference volume for entire chunk
        [RS1, CS1] = DetermineXYShiftsFBS(chunk_reg0, blurFactor, keepFactor, ref1); %calculate DFT shift between each volume and target vol
        chunk_reg1 = ApplyXYShiftsFBS(chunk_reg0,RS1,CS1); %apply the shift (need it for subsequent steps)

        % 3D DFT registration
        ref2 = mean(chunk_reg1, 4); %defineReference(chunk_reg1, chunkSize, refType);
        shifts = zeros(chunkLength(c),3);
        if sbxInfo.Nplane > 1
            parfor j = 1:chunkLength(c)
                shifts(j,:) = dftregistration3D(fftn(ref2), fftn(chunk_reg1(:,:,:,j)), 2); % dftregistration3D(fftn(ref2), fftn(vol), 2);
            end
        end
        RS2 = shifts(:,1)*regParam.binXY;
        CS2 = shifts(:,2)*regParam.binXY;
        ZS1 = shifts(:,3);

        % combine Row and Column Shifts from DFT registrations above
        RS(:,c) = {RS0*regParam.binXY + RS1*regParam.binXY + repmat(RS2',[sbxInfo.Nplane,1])};
        CS(:,c) = {CS0*regParam.binXY + CS1*regParam.binXY + repmat(CS2',[sbxInfo.Nplane,1])};
        ZS(:,c) = {ZS1'};

        % save reference files for stitching later
        ref_all{c} = ref2;
        RS0_all{c} = RS0;
        CS0_all{c} = CS0;
        RS1_all{c} = RS1;
        CS1_all{c} = CS1;
        RS2_all{c} = RS2;
        CS2_all{c} = CS2;

        waitbar(c/Nchunk, w); % H.iterate(1); %
    end
    close(w); % close(H); %
    toc

    % What are these used for?
    %{
    intermediate_shifts.RS0_all = [RS0_all{:}];
    intermediate_shifts.RS1_all = [RS1_all{:}];
    intermediate_shifts.RS2_all = vertcat(RS2_all{:}); % [RS2_all{:}];
    intermediate_shifts.CS0_all = [CS0_all{:}];
    intermediate_shifts.CS1_all = [CS1_all{:}];
    intermediate_shifts.CS2_all = vertcat(CS2_all{:}); %[CS2_all{:}];
    %}

    % Fix intra-chunk discontinuities
    if anchor == 0, anchor = round(Nchunk/2); end
    ref_final = ref_all{anchor};
    interchunk_shifts = zeros(Nchunk,3);
    if sbxInfo.Nplane > 1
        for j = 1:Nchunk
            interchunk_shifts(j,:) = dftregistration3D(fftn(ref_final),fftn(ref_all{j}),2);
        end
    end
    RS_chunk = interchunk_shifts(:,1)*regParam.binXY;
    CS_chunk = interchunk_shifts(:,2)*regParam.binXY;
    ZS_chunk = interchunk_shifts(:,3);
    RS = [RS{:}];  CS = [CS{:}];  ZS = [ZS{:}]; %convert local shift correction cells to matrix form

    %stretch the intra-chunk corrections to apply to every frame
    RS_chunk = imresize(RS_chunk',size(RS),'nearest');
    CS_chunk = imresize(CS_chunk',size(CS),'nearest');
    ZS_chunk = imresize(ZS_chunk',size(ZS),'nearest');

    % save the result transformations, and other relevant info
    fprintf('\nSaving %s', shiftPath)
    save(shiftPath, 'expt','regParam','chunkSize','anchor','keepFactor','blurFactor','RS','CS','ZS','RS_chunk','CS_chunk','ZS_chunk','ref_all', '-mat', '-v7.3'); % , 'intermediate_shifts'
else
    fprintf('\nLoading %s', shiftPath)
    load(shiftPath);
end

% SECOND STAGE - FORMERLY part of MakeSbxDFT
if ~exist('CS_final','var')
    % Combine results from first stage
    ZS_total = ZS+ZS_chunk;
    ZS_final = ZS_total - median(ZS_total);
    RS_total = RS+RS_chunk;
    RS_total = RS_total - median(RS_total);
    CS_total = CS+CS_chunk;
    CS_total = CS_total - median(CS_total);

    %work one volume at a time, do the zproj registration
    R_zproj = []; C_zproj = []; % R_final = []; C_final = [];
    if expt.Nplane > 1 %projToggle
        w = waitbar(0, 'applying initial shifts');
        zproj_raw = zeros(1, Nx_crop, Ny_crop, sbxInfo.totScan);
        tic
        for scan = 1:sbxInfo.totScan
            %read the volume
            raw_vol = readSBX(sbxInputPath, sbxInfo, scan, 1, refPMT, []); % usePMT
            %if sbxInfo.nchan == 1, raw_vol = reshape(raw_vol, sbxInfo.nchan, Nx, Ny, sbxInfo.Nplane); end
            %crop it based on edges
            raw_vol = raw_vol(regParam.edges(3)+1:end-regParam.edges(4),regParam.edges(1)+1:end-regParam.edges(2),:);
            %first do individual plane XY shifts
            reg_vol = nan(size(raw_vol));
            for z = flip(1:sbxInfo.Nplane)
                reg_vol(:,:,z) = imtranslate(raw_vol(:,:,z),[CS_total(z,scan), RS_total(z,scan)]);
            end
            if sbxInfo.Nplane > 1,  reg_vol = imtranslate(reg_vol, [0,0,ZS_final(scan)]);  end  % then do z shift
            zproj_raw(1,:,:,scan) = mean(reg_vol(:,:,zProjRange),3);
            waitbar(scan/sbxInfo.totScan);
        end
        delete(w);

        fprintf('\nRegistering z-projection... '); tic
        [zproj_mean, R_zproj, C_zproj] = zproj_reg(1, sbxInfo.totScan, usePMT, zProjRange, 'refchan',1, 'zproj_raw',zproj_raw); % refPMT
        if ~isempty(zprojPath)
            write2chanTiff(uint16(zproj_mean), zprojPath);  
        end
        R_zproj = transpose(repmat(R_zproj,1,sbxInfo.Nplane));
        C_zproj = transpose(repmat(C_zproj,1,sbxInfo.Nplane));
        RS_final = RS_total + R_zproj;
        CS_final = CS_total + C_zproj;
        toc
    else
        %fprintf('\nz-projection registration disabled\n');
        RS_final = RS_total;
        CS_final = CS_total;
    end
    fprintf('\nUpdating %s', shiftPath)
    save( shiftPath, 'R_zproj','C_zproj','CS_final','RS_final', 'ZS_final', '-append' ) % save the projection-based results into the existing dft shifts file  ,'zRange'
end

% Apply the translations one scan at a time, and write those results to the sbxdft file
if ~exist(expt.sbx.dft, 'file')
    if expt.Nplane == 1 && sbxInfo.nchan == 1
        data_type = 1;  % single-plane, single-color
    elseif expt.Nplane == 1 && sbxInfo.nchan > 1
        data_type = 2;  % single-plane, multi-color
    elseif expt.Nplane > 1 && sbxInfo.nchan == 1
        data_type = 3;  % multi-plane, single-color
    elseif expt.Nplane > 1 && sbxInfo.nchan > 1
        data_type = 4;  % multi-plane, multi-color
    end
    if Npmt > 1, usePMT = -1; end
    tic
    w = waitbar(0,'writing .sbxdft');
    rw = SbxWriter(expt.sbx.dft, sbxInfo, '.sbxdft', true);
    fprintf('\n     Writing corrected sbx file');
    tic
    for scan = 1:sbxInfo.totScan
        data_scan = readSBX(sbxInputPath, sbxInfo, scan, 1, usePMT, []); % load each scan
        % Reshape the scan data to apply imtranslate, and then reshape again to write the new sbxdft file
        switch data_type
            case 1 % single-plane, single-color (check) 
                data_scan = imtranslate(data_scan, [CS_final(1,scan), RS_final(1,scan)]);
                %data_scan = permute(data_scan, [4,1,2,3]); 
            case 2 % Single plane, multi-color
                data_scan = permute(data_scan, [2,3,4,1]); % [x,y,t,c]
                for chan = 1:2
                    data_scan(:,:,1,chan) = imtranslate(data_scan(:,:,1,chan), [CS_final(1,scan), RS_final(1,scan)]);
                end
                data_scan = permute(data_scan, [4,1,2,3]);
            case 3 % multi-plane, single-color
                for z = 1:sbxInfo.Nplane
                    data_scan(:,:,z) = imtranslate(data_scan(:,:,z), [CS_final(z,scan), RS_final(z,scan)]);
                end
            case 4
                % multi plane, multi color (check)
                data_scan = permute(data_scan, [2,3,4,5,1]); % [x,y,z,t,c]
                for chan = 1:2
                    for z = 1:sbxInfo.Nplane
                        data_scan(:,:,z,1,chan) = imtranslate(data_scan(:,:,z,1,chan), [CS_final(z,scan), RS_final(z,scan)]);
                    end
                end
                data_scan = permute(data_scan, [5,1,2,3,4]); % [c,x,y,z,t]
                data_scan = reshape(data_scan, [size(data_scan,[1,2,3]), prod(size(data_scan,[4,5]))]);
        end
        rw.write(data_scan);
        waitbar(scan/sbxInfo.totScan, w);
    end
    rw.delete;
    delete(w);
end
end