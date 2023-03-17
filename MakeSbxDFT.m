function zproj_mean = MakeSbxDFT(sbxPath, sbxInfo, shiftPath, varargin)
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'sbxInfo', @isstruct )
addRequired(IP, 'shiftPath', @ischar )
addOptional(IP, 'refChan', 'green', @ischar ) % for scanbox, PMT1 = green, 2 = red. -1 = both
addParameter(IP, 'edges',[0,0,0,0], @isnumeric); % [left, right, top, bottom]
addParameter(IP, 'zprojPath', '', @ischar)
addParameter(IP, 'proj', true, @islogical);
addParameter(IP, 'range', [0.25, 0.85], @isnumeric);
parse(IP, sbxPath, sbxInfo, shiftPath, varargin{:}); % tforms_optotune, 
edges = IP.Results.edges;
refChan = IP.Results.refChan;
[refPMT, ~] = DeterminePMT(refChan, sbxInfo); % PMT1 = green, PMT2 = red
projToggle = IP.Results.proj;
zprojPath = IP.Results.zprojPath;
proj_range = IP.Results.range;

[usePMT, ~] = DeterminePMT('both', sbxInfo);
if numel(usePMT) > 1, usePMT = -1; end


load(shiftPath,'-mat');
Nrow = 512; %sbxInfo.sz(1);  
Ncol = 796; %sbxInfo.sz(2); 
Nscan = size(CS,2);
if ~exist('CS_final','var')
    if sbxInfo.Nplane > 1
        zProjRange = round(proj_range(1)*sbxInfo.Nplane):round(proj_range(2)*sbxInfo.Nplane);
    else
        zProjRange = [1, 1];
    end
    if isempty( zprojPath ) 
        [fdir, fname ] = fileparts( sbxPath );
        zprojPath = sprintf('%s\\%s_zproj.tif', fdir, fname);
    end

    % Load results from CorrectData3D
    ZS_total = ZS+ZS_chunk;
    ZS_final = ZS_total - median(ZS_total);
    RS_total = RS+RS_chunk;
    RS_total = RS_total - median(RS_total);
    CS_total = CS+CS_chunk;
    CS_total = CS_total - median(CS_total);
    NrowCrop = Nrow - edges(3) - edges(4);
    NcolCrop = Ncol - edges(1) - edges(2);
    
    %work one volume at a time, do the zproj registration
    tic   
    R_zproj = []; C_zproj = []; % R_final = []; C_final = [];
    if projToggle
        zproj_raw = zeros(sbxInfo.nchan, NrowCrop, NcolCrop, Nscan);
        w = waitbar(0, 'applying initial shifts');
        for s = 1:Nscan
            %read the volume
            raw_vol = readSBX(sbxPath, sbxInfo, s, 1, usePMT, []); %pipe.imread(sbxPath, sbxInfo.Nplane*(i-1)+1, sbxInfo.Nplane, pmt,[]); %readSBX(path, info, k, N, pmt, optolevel)
            raw_vol = reshape(raw_vol, sbxInfo.nchan, Nrow, Ncol, sbxInfo.Nplane);
            %crop it based on edges
            raw_vol = raw_vol(:,edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:);      
            tempVol = cell(1,sbxInfo.nchan); reg_vol = zeros(size(raw_vol)); 
            for c = 1:sbxInfo.nchan
                %first do individual plane XY shifts 
                A = cell(sbxInfo.Nplane,1);
                parfor z = 1:sbxInfo.Nplane
                    A{z} = imtranslate(squeeze(raw_vol(c,:,:,z)),[CS_total(z,s), RS_total(z,s)]);
                end
                tempVol{c} = cat(3, A{:});

                if sbxInfo.Nplane > 1
                    tempVol{c} = imtranslate(tempVol{c}, [0,0,ZS_final(s)]); % then do z shift
                end
                reg_vol(c,:,:,:) = tempVol{c};
            end
            zproj_raw(:,:,:,s) = mean(reg_vol(:,:,:,zProjRange),4);

            waitbar(s/Nscan);
        end
        delete(w);

        fprintf('\nRegistering z-projection... '); tic
        [zproj_mean, R_zproj, C_zproj] = zproj_reg(1, Nscan, usePMT, zProjRange, 'refchan',refPMT, 'zproj_raw',zproj_raw);
        if ~isempty(zprojPath),  write2chanTiff(uint16(zproj_mean), zprojPath);  end
        R_zproj = transpose(repmat(R_zproj,1,sbxInfo.Nplane));
        C_zproj = transpose(repmat(C_zproj,1,sbxInfo.Nplane));
        RS_final = RS_total + R_zproj;
        CS_final = CS_total + C_zproj;
        toc
    else
        fprintf('\nz-projection registration disabled\n');
        RS_final = RS_total; 
        CS_final = CS_total;
    end
    fprintf('\nUpdating %s', shiftPath)
    save( shiftPath, 'R_zproj','C_zproj','CS_final','RS_final', 'ZS_final', '-append' ) % save the projection-based results into the existing dft shifts file  ,'zRange'
end

%{
subplot(3,1,1); 
imagesc( CS_final ); % reshape( normalize(CS_final(:)), size(CS_final)  )
title('Column shift');
caxis([-3,3]);

subplot(3,1,2); 
imagesc( RS_final ); % reshape( normalize(CS_final(:)), size(CS_final)  )
title('Row shift');
caxis([-3,3]);

subplot(3,1,3); 
imagesc( ZS_final ); % reshape( normalize(CS_final(:)), size(CS_final)  )
title('Z shift');
caxis([-3,3]);
impixelinfo
%}

% Load the sbxopt data
fprintf('\nLoading %s', sbxPath); tic;
sbx_data = readSBX(sbxPath, sbxInfo, 1, Nscan, usePMT, []); toc
% Apply 3D DFT translations
w = waitbar(0,'Applying 3D DFT translations'); tic;
if sbxInfo.nchan == 1
    for s = 1:Nscan
        for z = 1:sbxInfo.Nplane
            sbx_data(:,:,z,s) = imtranslate(sbx_data(:,:,z,s), [CS_final(z,s), RS_final(z,s)]); %first do XY shifts
        end
        sbx_data(:,:,:,s) = imtranslate(sbx_data(:,:,:,s), [0,0,ZS_final(s)]); % then do z shift
        waitbar(s/Nscan, w); 
    end
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2]), prod(size(sbx_data,[3,4]))] );    
else
    sbx_data = permute(sbx_data, [2,3,4,5,1]); % imtranslate expects spatial dimensions first
    for s = 1:Nscan
        for z = 1:sbxInfo.Nplane
            for c = 1:2
                sbx_data(:,:,z,s,c) = imtranslate(sbx_data(:,:,z,s,c), [CS_final(z,s), RS_final(z,s)]); %first do XY shifts
            end
        end
        sbx_data(:,:,:,s,c) = imtranslate(sbx_data(:,:,:,s,c), [0,0,ZS_final(s)]); % then do z shift
        waitbar(s/Nscan, w); 
    end
    sbx_data = permute(sbx_data, [5,1,2,3,4]); % regwriter expects color dim first
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2,3]), prod(size(sbx_data,[4,5]))] );    
end
delete(w);
toc

% Make the sbxdft file
fprintf('\nWriting sbxdft'); tic;  % , sbxPath
rw = SbxWriter(sbxPath, sbxInfo, '.sbxdft'); %pipe.io.RegWriter(sbxPath, sbxInfo, '.sbxdft'); %open regwriter object for writing sbxdft
rw.write(sbx_data);
rw.delete;
toc
end