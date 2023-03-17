function [volData, binLims, volPaths] = WriteSbxVolumeTif(sbxPath, sbxInfo, varargin)
checkInfo = @(x)(isstruct(x) || isempty(x));
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', checkInfo)
if isempty(sbxInfo)
    [pathDir, pathName, ~] = fileparts(sbxPath); % pathExt
    [~,infoPath] = FileFinder(pathDir, 'type','mat', 'criteria',@(x)(strcmp(x,pathName))); %
    fprintf('\nNo info structure was provided. Loading %s', infoPath{1})
    sbxInfo = MakeInfoStruct( infoPath{1} );
end
addParameter( IP, 'z', 1:sbxInfo.Nplane, @isnumeric )
addParameter( IP, 'chan', 'both', @ischar ) % 'green', 'red', 'both'. for scanbox, 1 = green, 2 = red. -1 = both
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'Nscan', -1, @isnumeric )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric ) % [left, right, top, bottom]
addParameter( IP, 'scale', 1, @isnumeric ) %
addParameter( IP, 'binT', 1, @isnumeric )
%addParameter( IP, 'zeros', true, @islogical )
addParameter( IP, 'dir', '', @ischar )
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'type', '', @ischar )
addParameter( IP, 'monochrome', true, @islogical ) % should we write bioformat tifs?
addParameter( IP, 'RGB', false, @islogical ) % should we write bioformat tifs?
addParameter( IP, 'verbose', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxPath, sbxInfo, varargin{:} );
% Image-related parameters
zGet = IP.Results.z;
writeChan = IP.Results.chan;  %1; %1 = red, 2 = green
firstScan = IP.Results.firstScan;
Nscan = IP.Results.Nscan;
edges = IP.Results.edges;
scaleFactor = IP.Results.scale;
binT = IP.Results.binT;
%allowZeros = IP.Results.zeros;
verbose = IP.Results.verbose;
outputToggle = nargout > 0;

% Determine output file path and name
monochrome = IP.Results.monochrome;
RGBtoggle = IP.Results.RGB;
overwrite = IP.Results.overwrite;
saveDir = IP.Results.dir;
saveName = IP.Results.name;
if isempty(saveName), [~,saveName] = fileparts(sbxPath); end
sbxType = IP.Results.type;
if isempty(saveDir), saveDir = strcat(fileparts(sbxPath), '\'); end
mkdir(saveDir);
if ~isempty(sbxType)
    nameRoot = sprintf('%s_%s', saveName, sbxType);
else
    nameRoot = sprintf('%s', saveName);
end
binLims = []; volData = [];
if outputToggle || monochrome || RGBtoggle
    % DETERMINE WHICH CHANNELS TO WRITE
    [usePMT, ~] = DeterminePMT(writeChan, sbxInfo);
    Npmt = numel(usePMT);
    if Npmt == 1
        pmtArg = usePMT;
        Ndim = 4;
    else
        pmtArg = -1;
        Ndim = 5;
    end

    % DETERMINE # OF SCANS, BINNING, ETC
    if Nscan < 0,  Nscan = sbxInfo.totScan;  end
    [binLims, Nbin, binLength] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',binT );

    % GET AND PROCESS IMAGING DATA
    if verbose, tic; end
    if binT == 1
        volData = readSBX(sbxPath, sbxInfo, firstScan, Nscan, pmtArg);
    else
        if verbose, fprintf('\nLoading %s (first scan = %i, Nscan = %i, pmt = %i, averaging every %i frames)    ', sbxPath, firstScan, Nscan, pmtArg, binT ); end
        volData = cell(1,Nbin);
        w = parfor_progressbar(Nbin, sprintf('Binning %s', sbxPath )); %waitbar(0, sprintf('Binning %s', sbxPath ));  
        tic;
        parfor b = 1:Nbin %flip(1:Nbin) % parfor is faster
            volChunk = readSBX(sbxPath, sbxInfo, binLims(b,1), binLength(b), pmtArg);
            volData{b} = mean(volChunk, Ndim);
            w.iterate(1); %waitbar((Nbin-b)/Nbin, w);
        end
        toc
        close(w); % delete(w);
        clearvars volChunk;
        volData = cellfun(@uint16, volData, 'UniformOutput',false);
        toc

        tic
        volData = cat(Ndim, volData{:});
        toc
    end
    % Crop and rearrange the data
    tic
    if Ndim == 4
        volData = volData(edges(3)+1:end-edges(4), edges(1)+1:end-edges(2), zGet, :, :);
        volData = permute(volData, [1,2,3,5,4]); % [x,y,z,c,t] -> [x,y,z,c,t]
    else
        volData = volData(:,edges(3)+1:end-edges(4), edges(1)+1:end-edges(2), zGet, :);
        volData = permute(volData, [2,3,4,1,5]); % [c,x,y,z,t] -> [x,y,z,c,t]
        volData = volData(:,:,:,[2,1],:); % GR -> RG
    end
    toc
    
    % Rescale
    if scaleFactor ~= 1
        scaleDim = [size(volData,1)/scaleFactor, size(volData,2)/scaleFactor, size(volData,3)];
        for c = Npmt:-1:1
            for s = flip(1:size(volData,5))
                rescaleData(:,:,:,c,s) = imresize3( volData(:,:,:,c,s), scaleDim );
            end
        end
        volData = rescaleData;
        clearvars rescaleData;
    end

    % WRITE TIFFS
    chanName = ["red","green"];
    if monochrome
        chanTifPath = cell(1,2);
        for chan = 1:size(volData,4)
            chanTifPath{chan} = sprintf('%s%s_vol_%s.tif', saveDir, nameRoot, chanName{chan});
            if overwrite || ~exist(chanTifPath{chan},'file')
                fprintf('\nWriting Bio-Formats tiff %s', chanTifPath{chan}); tic;
                bfsave(uint16(volData(:,:,:,chan,:)), chanTifPath{chan}); toc
            end
            toc
        end
    end
    rgbPath = sprintf('%s%s_vol_RGB.tif', saveDir, nameRoot); 
    if RGBtoggle
        if overwrite || ~exist(rgbPath,'file')
            tic
            for chan = flip(1:size(volData,4))
                chanLower = prctile(reshape(volData(:,:,:,chan,:),[],1), 5);
                chanUpper = max(reshape(volData(:,:,:,chan,:),[],1)); %prctile(stackChan{chan}(:), 1);
                fprintf('\nRescaling %s channel : [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper);
                rgbData(:,:,:,chan,:) = uint8(rescale(volData(:,:,:,chan,:), 0, 2^8-1, 'inputMin',chanLower, 'inputMax',chanUpper)); toc
            end
            fprintf('\nWriting Bio-Formats tiff %s', rgbPath);
            bfsave(rgbData, rgbPath); % , 'BigTiff',true
            toc
        end
    end
    volPaths = [chanTifPath, {rgbPath}];
end
end