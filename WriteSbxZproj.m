function [projData, binLims, projPaths, rgbData] = WriteSbxZproj(sbxPath, sbxInfo, varargin) %  projPath,  firstScan, Nscan, pmt,
checkInfo = @(x)(isstruct(x) || isempty(x) || ischar(x));
checkz = @(z)(isnumeric(z) || iscell(z));
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar ) % path to the .sbx+ file to be projected
% Get the metadata associated with this sbx file
addRequired( IP, 'sbxInfo', checkInfo ) % @isstruct
if ischar(sbxInfo)
    sbxInfo = MakeInfoStruct( sbxInfo );
elseif isempty(sbxInfo)
    [pathDir, pathName, ~] = fileparts(sbxPath); % pathExt
    [~,infoPath] = FileFinder(pathDir, 'type','mat', 'criteria',@(x)(strcmp(x,pathName))); %
    if ~isempty(infoPath)
        fprintf('\nNo info structure was provided. Loading %s', infoPath{1})
        sbxInfo = MakeInfoStruct( infoPath{1} );
    else
        error('\nUnable to find sbxInfo file!')
    end
end
addParameter( IP, 'projType', 'mean', @ischar)
addParameter( IP, 'z', 1:sbxInfo.Nplane, checkz ) % @isnumeric
addParameter( IP, 'edges', [0, 0, 0, 0], @isnumeric ) % [left, right, top, bottom]
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'Nscan', -1, @isnumeric )
addParameter( IP, 'chan', 'green', @ischar )
addParameter( IP, 'scale', 1, @isnumeric ) %
addParameter( IP, 'binT', 1, @isnumeric )
addParameter( IP, 'write', true, @islogical )
addParameter( IP, 'dir', '', @ischar )
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'sbxType', '', @ischar )
addParameter( IP, 'monochrome', true, @islogical ) % write a monochrome movie for each color?
addParameter( IP, 'RGB', false, @islogical ) % write an RGB movie?
addParameter( IP, 'overwrite', false, @islogical ) % for scanbox, 1 = green, 2 = red. -1 = both
parse( IP, sbxPath, sbxInfo, varargin{:} ); % sbxInfo, projPath,  mouse, exptDate,
% Projection parameters
projType = IP.Results.projType;
zProj = IP.Results.z;
if isnumeric(zProj), zProj ={zProj};  end
Nz = numel(zProj);
edges = IP.Results.edges;
projChan = IP.Results.chan;
scaleFactor = IP.Results.scale;
firstScan = IP.Results.firstScan;
Nscan = IP.Results.Nscan;
if Nscan == -1, Nscan = sbxInfo.totScan; end % sbxInfo.Nscan
% set up temporal binning
binT = IP.Results.binT;
if binT ~= 1 && rem(Nscan, binT) ~= 0
    cutScans = Nscan-binT*floor(Nscan/binT);
    fprintf('\nNumber of scans (%i) is not divisible by %i: cutting off first %i scans',Nscan, binT, cutScans);
    firstScan = firstScan+cutScans;
end
[binLims, ~, ~] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',binT, 'allowPartial',false);
projLength = binLims(end)-binLims(1)+1;
% Determine output file paths and names
monochrome = IP.Results.monochrome;
RGB = IP.Results.RGB;
overwrite = IP.Results.overwrite;
saveDir = IP.Results.dir;
if isempty(saveDir), saveDir = strcat(fileparts(sbxPath), '\'); end
mkdir(saveDir);
saveName = IP.Results.name;
if isempty(saveName), [~,saveName] = fileparts(sbxPath); end
sbxType = IP.Results.sbxType;

if nargout == 0 && ~(monochrome || RGB)
    fprintf('\nNo output was requested - skipping\n');
else
    % Determine which PMTs/channels to use
    chanName = {'red','green'}; %
    chanInd = [2,1]; %pmtInd = 1:2;
    [usePMT, ~] = DeterminePMT(projChan, sbxInfo);
    useChan = sort(chanInd(usePMT));
    Nchan = numel(useChan);
    if Nchan == 2 && sbxInfo.nchan == 2
        projDim = 4;
        projPMT = -1;
    else
        projDim = 3;
        projPMT = usePMT;
        RGB = false;
    end
    % Set up output file names, check if the files already exist, and load them if asked
    chanProjPath = cell(Nz,2); chanProjExists = false(Nz,2); chanData = cell(Nz,2); rgbPath = cell(Nz,1); projData = cell(Nz,1); rgbData = cell(Nz,1);
    fix(clock)
    for z = 1:Nz
        if ~isempty(sbxType)
            nameRoot = sprintf('%s_%s_%s_z%i-%i_', saveName, sbxType, projType, zProj{z}(1), zProj{z}(end));
        else
            nameRoot = sprintf('%s_%s_z%i-%i_', saveName, projType, zProj{z}(1), zProj{z}(end));
        end
        for chan = useChan
            chanProjPath{z,chan} = sprintf('%s%s%s.tif', saveDir, nameRoot, chanName{chan}); % chanName
            if exist(chanProjPath{z,chan},'file') && ~overwrite
                chanProjExists(z,chan) = true;
                fprintf('\nLoading %s', chanProjPath{z,chan}); %if verbose,  end
                chanData{z,chan} = tiffreadVolume(chanProjPath{z,chan}); % loadtiff(chanProjPath{z,chan});
            end
        end
        rgbPath{z} = sprintf('%s%sRGB.tif', saveDir, nameRoot);
        if exist(rgbPath{z}, 'file')
            fprintf('\nLoading %s', rgbPath{z});
            rgbData{z} = loadtiff(rgbPath{z});
        end
    end
    if RGB
        projPaths = [chanProjPath, rgbPath];
    else
        projPaths = chanProjPath(:,useChan);
    end

    % Load each scan and perform the projection(s)
    if ~all(chanProjExists(:,useChan), 'all') || overwrite
        tic
        k = 0;
        %vol_proj = [];
        for scan = binLims(1):binLims(end) %  par   firstScan+Nscan-1
            k = k+1;
            vol = readSBX(sbxPath, sbxInfo, scan, 1, projPMT, []); % load each volume scan
            for z = 1:Nz % par is slower
                if Nchan == 1
                    vol_proj = vol(:,:,zProj{z});
                else
                    vol_proj = vol(:,:,:,zProj{z});
                end
                vol_proj(vol_proj==0) = NaN; % suppress zeros from projection
                if strcmpi(projType,'mean')
                    projData{z}{k} = squeeze( mean(vol_proj,projDim,'omitnan') );
                elseif strcmpi(projType,'max')
                    projData{z}{k} = squeeze( max(vol_proj,[],projDim,'omitnan') );
                elseif strcmpi(projType,'median')
                    projData{z}{k} = squeeze( median(vol_proj,projDim,'omitnan') );
                else
                    disp('invalid projection type');
                end
                projData{z}{k}(isnan(projData{z}{k})) = 0;
            end
        end
        clearvars vol_proj vol
        % Post-process each projection and write tifs (optional)
        for z = 1:Nz
            % Concatenate projections
            projData{z} = cat(projDim, projData{z}{:});
            if ndims(projData{z}) == 4,  projData{z} = permute(projData{z}, [2,3,4,1]); end % permute(projData, [2,3,1,4])
            projData{z} = projData{z}(edges(3)+1:end-edges(4),  edges(1)+1:end-edges(2), :, :);% crop edges
            % Spatial averaging
            if scaleFactor ~= 1
                scaleSize = round([size(projData{z},1)/scaleFactor, size(projData{z},2)/scaleFactor, size(projData{z},3)]);
                if size(projData{z},4) > 1 %sbxInfo.nchan > 1
                    scaleData = cell(1,2);
                    for chan = useChan,  scaleData{chan} = imresize3( projData{z}(:,:,:,chan), scaleSize);  end
                    projData{z} = cat(4, scaleData{useChan});
                    clearvars scaleData;
                else
                    projData{z} = imresize3( projData{z}, scaleSize);
                end
            end
            % Temporal averaging
            if binT ~= 1
                projData{z} = reshape(projData{z}, size(projData{z},1), size(projData{z},2), binT, [], Nchan);
                projData{z} = squeeze(mean(projData{z}, 3));
            end

            % Switch channel order GRB -> RGB
            if Nchan > 1, projData{z} = flip(projData{z}, 4); end
            % write the monochrome projections to tif
            if monochrome
                if Nchan > 1
                    for chan = useChan
                        fprintf('\nWriting %s', chanProjPath{z,chan})
                        WriteTiff(uint16(projData{z}(:,:,:,chan)), chanProjPath{z,chan});
                    end
                else
                    fprintf('\nWriting %s', chanProjPath{z,useChan})
                    WriteTiff(uint16(projData{z}), chanProjPath{z,useChan});
                end
            end
            % Write or load the RGB projections (reduced to 8bit chan)
            if RGB && Nchan > 1
                redChan = projData{z}(:,:,:,1);
                redChan = uint8(rescale(redChan, 0, 255, 'InputMin',prctile(redChan(:),1), 'InputMax',prctile(redChan(:),99)));
                greenChan = projData{z}(:,:,:,2);
                greenChan = uint8(rescale(greenChan, 0, 255, 'InputMin',prctile(greenChan(:),1), 'InputMax',prctile(greenChan(:),99)));
                rgbData{z} = cat(4, redChan, greenChan);
                clearvars redChan greenChan;
                rgbData{z}(end,end,end,3) = 0;  % need a blue channel for RGB data
                WriteTiff(rgbData{z}, rgbPath{z});
            end
            toc
        end
    else
        % Combine loaded monochrome data into projData
        for z = 1:Nz
            projData{z} = cat(projDim, chanData{z,:});
        end
    end
    % Convert cell to array, if only one zProj was given
    if Nz == 1
        projData = projData{1};
    end
end