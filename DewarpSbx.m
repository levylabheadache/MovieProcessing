function tforms_optotune = DewarpSbx( expt, sbxInfo, regParam, varargin)
IP = inputParser;
addRequired(IP, 'expt', @isstruct )
addRequired(IP, 'sbxInfo', @isstruct )
addRequired(IP, 'regParam', @isstruct )
addParameter(IP, 'type', 'rigid', @ischar) %options are 'affine', 'rigid', or 'none' 'affine'
addParameter(IP, 'overwrite', false, @islogical)
addParameter(IP, 'show', false, @islogical)
parse(IP, expt, sbxInfo, regParam, varargin{:});
regType = IP.Results.type;
overwrite = IP.Results.overwrite;
show = IP.Results.show;
optPath = expt.sbx.cat; optPath(end-2:end) = 'opt';

if ~exist(optPath,'file') || overwrite
    [refPMT, ~] = DeterminePMT(regParam.refChan, sbxInfo); % PMT1 = green, PMT2 = red  refPMTname
    tformPath = sprintf('%s%s_dewarp.mat', expt.dir, expt.name );
    % CALCULATE OPTOTUNE WARPING
    if ~exist(tformPath, 'file') || overwrite
        % Load reference volume, and crop edges
        cropRows = regParam.edges(3)+1:sbxInfo.sz(1)-regParam.edges(4); Nrow = numel(cropRows);
        cropCol = regParam.edges(1)+1:sbxInfo.sz(2)-regParam.edges(2); Ncol = numel(cropCol);
        raw_ref = readSBX(expt.sbx.cat, sbxInfo, regParam.refScan(1), numel(regParam.refScan), refPMT, []);
        raw_ref = raw_ref(cropRows, cropCol, :);
        raw_ref = reshape(raw_ref, Nrow, Ncol, sbxInfo.Nplane, []);

        % Define a reference volume
        ref_vol = squeeze( mean(raw_ref,4)  ); % squeeze(median(raw_ref,4));

        %calculate optotune warping transformation from the mean volume
        fdir = fileparts(expt.sbx.cat);
        if strcmp(regType,'affine')
            tforms_optotune = MultiStackReg_Fiji_affine(ref_vol, fdir, sbxInfo.Nplane); % _2
            for z = 1:sbxInfo.Nplane, tforms_optotune(z).T([2,4]) = 0;  end % suppress shearing
        elseif strcmp(regType, 'rigid')
            tforms_optotune = OptoAlign_rigid(ref_vol); % align_vol
        elseif strcmpi(regType, 'none')
            tforms_optotune = repmat(affine2d(eye(3)),[1,sbxInfo.Nplane]);
        else
            fprintf('Invalid registration type for optotune correction');
            return;
        end

        %save the transformations for later
        fprintf('\nSaving %s', tformPath)
        save(tformPath,'tforms_optotune','expt', 'sbxInfo', 'regParam','regType'); %strcat(fdir,filesep,'tforms_optotune.mat')
    else
        fprintf('\nLoading %s', tformPath)
        load(tformPath,'tforms_optotune');
    end
    
    % Show the results (optional)
    if show
        for z = 1:sbxInfo.Nplane
            deformStruct.trans_x(z) = tforms_optotune(z).T(3,1);
            deformStruct.trans_y(z) = tforms_optotune(z).T(3,2);
            deformStruct.scale_x(z) = tforms_optotune(z).T(1,1);
            deformStruct.scale_y(z) = tforms_optotune(z).T(2,2);
            deformStruct.shear_x(z) = tforms_optotune(z).T(1,2);
            deformStruct.shear_y(z) = tforms_optotune(z).T(2,1);
        end
        opt = {[0.05,0.07], [0.07,0.04], [0.2,0.2]};  % {[vert, horz], [bottom, top], [left, right] }
        %close all;  clearvars sp;
        figure('WindowState','maximized', 'Color','w', 'PaperOrientation','landscape');
        sp(1) = subtightplot(6,1,1,opt{:});  plot( deformStruct.trans_x ); ylabel('X Trans'); title( sprintf('%s registration, ref size = %i scans', regType, refSize) ); %scale = %i, , scaleFactor
        sp(2) = subtightplot(6,1,2,opt{:});  plot( deformStruct.trans_y ); ylabel('Y Trans');
        sp(3) = subtightplot(6,1,3,opt{:});  plot( deformStruct.scale_x ); ylabel('X Scale');
        sp(4) = subtightplot(6,1,4,opt{:});  plot( deformStruct.scale_y ); ylabel('Y Scale');
        sp(5) = subtightplot(6,1,5,opt{:});  plot( deformStruct.shear_x ); ylabel('X Shear');
        sp(6) = subtightplot(6,1,6,opt{:});  plot( deformStruct.shear_y ); ylabel('Y Shear'); xlabel('Optotune Plane');
        linkaxes(sp,'x');
    end

    % Apply the transforms and generate sbxopt
    [chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.totScan, 'N',10);
    w = waitbar(0,'Applying optotune corrections');
    rw = SbxWriter(optPath, sbxInfo, '.sbxopt', true);
    imRef = imref2d([sbxInfo.sz(1), sbxInfo.sz(2)]);
    tic
    if sbxInfo.nchan == 1 
        for chunk = 1:Nchunk
            data_chunk = readSBX(expt.sbx.cat, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), refPMT, []); % [rows, col, planes, scans]
            for scan = 1:chunkLength(chunk)
                for z = 1:sbxInfo.Nplane
                    data_chunk(:,:,z,scan) = imwarp(data_chunk(:,:,z,scan), tforms_optotune(z), 'OutputView',imRef);
                end
            end
            data_chunk = reshape(data_chunk, [size(data_chunk,[1,2]), prod(size(data_chunk,[3,4]))]);
            rw.write(data_chunk);
            waitbar(chunk/Nchunk, w);
        end
    else
        for chunk = 1:Nchunk
            data_chunk = readSBX(expt.sbx.cat, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), -1, []);
            data_chunk = permute(data_chunk, [2,3,4,5,1]);
            for scan = 1:chunkLength(chunk)
                for chan = 1:2
                    for z = 1:sbxInfo.Nplane
                        data_chunk(:,:,z,scan,chan) = imwarp(data_chunk(:,:,z,scan,chan), tforms_optotune(z), 'OutputView',imRef);
                    end
                end
            end
            data_chunk = permute(data_chunk, [5,1,2,3,4]);
            data_chunk = reshape(data_chunk, [size(data_chunk,[1,2,3]), prod(size(data_chunk,[4,5]))]);
            rw.write(data_chunk);
            waitbar(chunk/Nchunk, w);
        end
    end
    toc
    rw.delete;
    delete(w);
end
end