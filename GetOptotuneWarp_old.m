function tforms_optotune = GetOptotuneWarp( expt, sbxPath, sbxInfo, varargin)
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'sbxInfo', @isstruct )
addParameter( IP, 'chan', 'green', @ischar ) % 'green', 'red', 'both'. for scanbox, 1 = green, 2 = red. -1 = both
addParameter(IP, 'firstRefScan', 3, @isnumeric);
addParameter(IP, 'Nref', 30, @isnumeric); % # of scans to use for reference volume
addParameter(IP, 'edges',[0,0,0,0], @isnumeric); % [left, right, top, bottom]
addParameter(IP, 'type', 'rigid', @ischar); %options are 'affine', 'rigid', or 'none' 'affine'
addParameter(IP, 'dir', '', @ischar);
%addParameter(IP, 'save', true, @islogical); % directory to save to
addParameter(IP, 'show', false, @islogical);
parse(IP, sbxPath, sbxInfo, varargin{:});
refChan = IP.Results.chan;
firstRefScan = IP.Results.firstRefScan;
refSize = IP.Results.Nref;
edges = IP.Results.edges;
regType = IP.Results.type;
fovDir = IP.Results.dir;
%saveToggle = IP.Results.save;
%saveToggle = ~isempty(saveDir);
show = IP.Results.show;

%sbxPathParts = strsplit(sbxPath, '\');
%fovDir = [strjoin(sbxPathParts(1:end-2), '\'),'\'];   
tformPath = sprintf('%s%stforms_optotune.mat',fovDir, ); %'D:\2photon\CGRPAi6-1\220816_FOV1\tforms_optotune.mat'; % strcat(fdir,filesep,'tforms_optotune.mat');
if ~exist(tformPath, 'file')
    % CALCULATE OPTOTUNE WARPING
    cropRows = edges(3)+1:sbxInfo.sz(1)-edges(4); Nrow = numel(cropRows);
    cropCol = edges(1)+1:sbxInfo.sz(2)-edges(2); Ncol = numel(cropCol);

    % Load n volumes, and crops the edges - NOTE: very first frame of the sbx is sometimes bad, better to start from the second scan
    if strcmpi(refChan, 'green')
        raw_ref = readSBX(sbxPath, sbxInfo, firstRefScan, refSize, 1, []);
    elseif strcmpi(refChan, 'red')
        raw_ref = readSBX(sbxPath, sbxInfo, firstRefScan, refSize, 2, []);
    end
    raw_ref = raw_ref(cropRows, cropCol, :);
    raw_ref = reshape(raw_ref, Nrow, Ncol, sbxInfo.Nplane, []);

    % Define a reference volume
    ref_vol = squeeze( mean(raw_ref,4)  ); % squeeze(median(raw_ref,4));

    %calculate optotune warping transformation from the mean volume
    fdir = fileparts(sbxPath);
    if strcmp(regType,'affine')
        tforms_optotune = MultiStackReg_Fiji_affine(ref_vol,fdir,sbxInfo.Nplane); % _2
        for z = 1:sbxInfo.Nplane, tforms_optotune(z).T([2,4]) = 0;  end % suppress shearing
    elseif strcmp(regType, 'rigid')
        tforms_optotune = OptoAlign_rigid(ref_vol); % align_vol
        %{
        tforms_optotune = MultiStackReg_Fiji_rigid(ref_vol,fdir,sbxInfo.Nplane); % rigid still returns an affine transformation matrix?
        for z = 1:sbxInfo.Nplane
            tforms_optotune(z).T([2,4,7,8]) = 0; % suppress shearing
            tforms_optotune(z).T([1,5,9]) = 1; % suppress scaling
        end
        %}
    elseif strcmpi(regType, 'none')
        tforms_optotune = repmat(affine2d(eye(3)),[1,sbxInfo.Nplane]);
    else
        fprintf('Invalid registration type for optotune correction');
        return;
    end

    %save the transformations for later
    if saveToggle
        fprintf('\nSaving %s', tformPath)
        save(tformPath,'tforms_optotune','firstRefScan','refSize','edges','refChan','regType'); %strcat(fdir,filesep,'tforms_optotune.mat')
    end
else
    fprintf('\nLoading %s', tformPath)
    load(tformPath,'tforms_optotune');
end
    
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

% Write optotune-corrected sbx file (sbxopt)
fprintf('\nLoading %s', sbxPath); tic;
sbx_data = readSBX(sbxPath, sbxInfo, 1, sbxInfo.Nscan, -1, []); toc; % load the raw data
w = waitbar(0,'Applying optotune corrections'); tic;
fprintf('\n Applying optotune corrections'); 
if sbxInfo.nchan == 1
    imRef = imref2d(size(sbx_data,[1,2]));
    for s = 1:sbxInfo.Nscan
        for z = flip(1:sbxInfo.Nplane)
            sbx_data(:,:,z,s) = imwarp( sbx_data(:,:,z,s), tforms_optotune(z), 'OutputView',imRef); % apply optotune dewarping
        end
        waitbar(s/sbxInfo.Nscan, w); 
    end
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2]), prod(size(sbx_data,[3,4]))]);
else
    imRef = imref2d(size(sbx_data,[2,3]));
    sbx_data = permute(sbx_data, [2,3,4,5,1]); % imtranslate expects spatial dimensions first
    % Apply the transformations to each color, scan by scan
    for s = 1:sbxInfo.Nscan
        for c = 1:2
            for z = flip(1:sbxInfo.Nplane) % 
                sbx_data(:,:,z,s,c) = imwarp( sbx_data(:,:,z,s,c), tforms_optotune(z), 'OutputView',imRef); 
            end
        end
        waitbar(s/sbxInfo.Nscan, w); 
    end
    sbx_data = permute(sbx_data, [5,1,2,3,4]); % regwriter expects color dim first
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2,3]), prod(size(sbx_data,[4,5]))] );    
end
delete(w);
% Make the sbxopt file
optPath = sbxPath; optPath(end-2:end) = 'opt';
fprintf('\nWriting %s', optPath); tic;
rw = SbxWriter(sbxPath, sbxInfo, '.sbxopt', true); % pipe.io.RegWriter(sbxPath, sbxInfo, '.sbxopt', true);
rw.write(sbx_data);
rw.delete;
toc
end