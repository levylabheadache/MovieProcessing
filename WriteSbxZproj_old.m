function [projData, binLims, projPaths] = WriteSbxZproj(sbxPath, sbxInfo, varargin) %  projPath,  firstScan, Nscan, pmt, 
checkInfo = @(x)(isstruct(x) || isempty(x));
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar ) % path to the .sbx+ file to be projected
addRequired( IP, 'sbxInfo', checkInfo ) % @isstruct
if isempty(sbxInfo)
    [pathDir, pathName, ~] = fileparts(sbxPath); % pathExt
    [~,infoPath] = FileFinder(pathDir, 'type','mat', 'criteria',@(x)(strcmp(x,pathName))); % 
    fprintf('\nNo info structure was provided. Loading %s', infoPath{1})
    sbxInfo = MakeInfoStruct( infoPath{1} );
end
addParameter( IP, 'projType', 'mean', @ischar)
addParameter( IP, 'z', 1:sbxInfo.Nplane, @isnumeric )
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
[binLims, ~, ~] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',binT, 'allowPartial',true);
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
if ~isempty(sbxType)
    nameRoot = sprintf('%s_%s_%s_z%i-%i_', saveName, sbxType, projType, zProj(1), zProj(end)); 
else
    nameRoot = sprintf('%s_%s_z%i-%i_', saveName, projType, zProj(1), zProj(end));
end
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
    end
    % Check if the data already exists, and load it if asked
    chanProjPath = cell(2,1); chanProjExists = false(2,1); chanData = cell(2,1);
    for chan = useChan
        chanProjPath{chan} = sprintf('%s%s%s.tif', saveDir, nameRoot, chanName{chan}); % chanName
        if exist(chanProjPath{chan},'file') && ~overwrite
            chanProjExists(chan) = true;
            fprintf('\nLoading %s', chanProjPath{chan}); %if verbose,  end
            chanData{chan} = loadtiff(chanProjPath{chan});
        end
    end
    rgbPath = sprintf('%s%sRGB.tif', saveDir, nameRoot);
    % Load each scan and perform the projection
    tic
    if ~all(chanProjExists(useChan)) || overwrite
        % load each volume scan (parallelized) and z-project
        A = cell(projLength,1); % firstScan+Nscan
        %h = parfor_progressbar(Nscan,'z-projecting...');
        k = 0; % tic
        for scan = binLims(1):binLims(end) % par   firstScan+Nscan-1
            k = k+1;
            vol = readSBX(sbxPath, sbxInfo, scan, 1, projPMT, []);
            if Nchan == 1
                vol = vol(:,:,zProj);
            else
                vol = vol(:,:,:,zProj);
            end
            vol(vol==0) = NaN;
            if strcmpi(projType,'mean')
                slice = squeeze( mean(vol,projDim,'omitnan') );
            elseif strcmpi(projType,'max')
                slice = squeeze( max(vol,[],projDim,'omitnan') );
            elseif strcmpi(projType,'median')
                slice = squeeze( median(vol,projDim,'omitnan') );
            else
                disp('invalid projection type');
            end
            slice(isnan(slice)) = 0;
            A{k} = slice;
            %h.iterate(1);
        end
        %delete(h); %toc
        % concatenate each projected volume
        projData = cat(projDim, A{:});
        clearvars A;
        if ndims(projData) == 4,  projData = permute(projData, [2,3,4,1]); end % permute(projData, [2,3,1,4])
        projData = projData(edges(3)+1:end-edges(4),  edges(1)+1:end-edges(2), :, :);% crop edges

        % Spatial averaging
        if scaleFactor ~= 1
            scaleSize = round([size(projData,1)/scaleFactor, size(projData,2)/scaleFactor, size(projData,3)]);
            if size(projData,4) > 1 %sbxInfo.nchan > 1
                scaleData = cell(1,2);
                for chan = useChan,  scaleData{chan} = imresize3( projData(:,:,:,chan), scaleSize);  end
                projData = cat(4, scaleData{useChan});
                clearvars scaleData;
            else
                projData = imresize3( projData, scaleSize);
            end
        end

        % Temporal averaging
        if binT ~= 1
            projData = reshape(projData, size(projData, 1), size(projData, 2), binT, [], Nchan);
            projData = squeeze(mean(projData, 3));
            %{
            if Nchan == 1
                projData = reshape(projData, size(projData, 1), size(projData, 2), binT, [], Nchan);
                projData = squeeze(mean(projData, 3));
            else
                projData = reshape(projData, size(projData, 1), size(projData, 2), binT, [], 2); %reshape(projData, size(projData, 1), size(projData, 2), [], binT, 2);
                projData = squeeze(mean(projData, 3));
            end
            %}
        end

        % switch chan order GRB -> RGB
        if Nchan > 1, projData = flip(projData, 4); end 
        toc
        % write the projection to a tif file
        if monochrome
            if Nchan > 1
                for chan = useChan
                    fprintf('\nWriting %s', chanProjPath{chan})
                    WriteTiff(uint16(projData(:,:,:,chan)), chanProjPath{chan});
                end
            else
                fprintf('\nWriting %s', chanProjPath{useChan})
                WriteTiff(uint16(projData), chanProjPath{useChan});
            end
        end
        if RGB && Nchan > 1
            redChan = projData(:,:,:,1);
            redChan = uint8(rescale(redChan, 0, 255, 'InputMin',prctile(redChan(:),1), 'InputMax',prctile(redChan(:),99)));
            greenChan = projData(:,:,:,1);
            greenChan = uint8(rescale(greenChan, 0, 255, 'InputMin',prctile(greenChan(:),1), 'InputMax',prctile(greenChan(:),99)));
            rgbData = cat(4, redChan, greenChan);
            clearvars redChan greenChan;
            rgbData(end,end,end,3) = 0;  % need a blue channel for RGB data
            WriteTiff(rgbData, rgbPath);
        end
    else
        projData = cat(4, chanData{useChan});
    end
    projPaths = [chanProjPath', {rgbPath}];
    toc
end
end