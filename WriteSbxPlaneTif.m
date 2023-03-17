function [stackOut, stackChan, binLims, projPaths] = WriteSbxPlaneTif(sbxPath, sbxInfo, z, varargin) 
checkInfo = @(x)(isstruct(x) || isempty(x) || ischar(x));
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', checkInfo ) 
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
addRequired( IP, 'z', @isnumeric)
addParameter( IP, 'chan', 'both', @ischar ) % 'green', 'red', 'both'. for scanbox, 1 = green, 2 = red. -1 = both
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'Nscan', -1, @isnumeric )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric ) % [left, right, top, bottom]
addParameter( IP, 'scale', 1, @isnumeric ) %
addParameter( IP, 'binT', 1, @isnumeric ) 
addParameter( IP, 'zeros', true, @islogical )
addParameter( IP, 'dir', '', @ischar ) 
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'sbxType', '', @ischar ) 
addParameter( IP, 'monochrome', true, @islogical ) % should we write individual monochrome tifs from each channel?
addParameter( IP, 'RGB', true, @islogical ) % should we write an RGB tif?
addParameter( IP, 'rescale', false, @islogical )
addParameter( IP, 'verbose', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxPath, sbxInfo, z, varargin{:} ); 
% Image-related parameters
writeChan = IP.Results.chan;  %1; %1 = red, 2 = green
firstScan = IP.Results.firstScan;
Nscan = IP.Results.Nscan;
edges = IP.Results.edges;
scaleFactor = IP.Results.scale;
binT = IP.Results.binT;
rescaleIntToggle = IP.Results.rescale;
allowZeros = IP.Results.zeros;
verbose = IP.Results.verbose;

outputToggle = nargout > 0;
if verbose, tic; end
% Determine output file path and name
monochrome = IP.Results.monochrome;
RGBtoggle = IP.Results.RGB;
overwrite = IP.Results.overwrite;
saveDir = IP.Results.dir;  
saveName = IP.Results.name;
if isempty(saveName), [~,saveName] = fileparts(sbxPath); end
sbxType = IP.Results.sbxType;
if isempty(saveDir), saveDir = strcat(fileparts(sbxPath), '\'); end
mkdir(saveDir);
if ~isempty(sbxType)
    nameRoot = sprintf('%s_%s', saveName, sbxType); 
else
    nameRoot = sprintf('%s', saveName);
end

% DETERMINE WHICH CHANNELS TO WRITE
pmtName = {'green','red'};
[usePMT, ~] = DeterminePMT(writeChan, sbxInfo);
Npmt = numel(usePMT);

% DETERMINE # OF SCANS, BINNING, ETC
if Nscan < 0,  Nscan = sbxInfo.totScan;  end
[binLims, Nbin, binLength] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',binT ); 

% GET DATA FROM EACH CHANNEL AND WRITE MONOCHROME TIFFS
stackChan = cell(1,2); pmtProjPath = cell(1,2);
for pmt = usePMT
    pmtProjPath{pmt} = sprintf('%s%s_Z%02d_%s.tif', saveDir, nameRoot, z, pmtName{pmt});
    fileExists = exist(pmtProjPath{pmt},'file');
    if outputToggle && ~overwrite && fileExists
        % Load stackChan from tif
        if verbose, fprintf('\nLoading %s', pmtProjPath{pmt});  end
        stackChan{pmt} = tiffreadVolume(pmtProjPath{pmt});  % loadtiff( pmtProjPath{pmt} );
        if ~allowZeros,  stackChan{pmt}(stackChan{pmt} == 0) = NaN;  end % Set zeros to NaN to avoid averaging over blank frames elsewhere
    elseif ~(~outputToggle && ~overwrite && fileExists)
        if binT == 1
            % Load, crop and resize the sbx data
            if verbose, fprintf('\nLoading %s (plane %i, first scan = %i, Nscan = %i, pmt = %i)    ', sbxPath, z, firstScan, Nscan, pmt ); end
            if sbxInfo.Nplane > 1
                stackChan{pmt} = readSBX(sbxPath, sbxInfo, firstScan, Nscan, pmt, z); 
            else
                stackChan{pmt} = readSBX(sbxPath, sbxInfo, firstScan, Nscan, pmt, []); 
            end
        else
            if verbose, fprintf('\nLoading %s (plane %i, first scan = %i, Nscan = %i, pmt = %i, averaging every %i frames)    ', sbxPath, z, firstScan, Nscan, pmt, binT ); end 
            w = waitbar(0, sprintf('Binning %s', sbxPath ));
            for c = flip(1:Nbin)
                if sbxInfo.Nplane > 1
                    tempChunk = readSBX(sbxPath, sbxInfo, binLims(c,1), binLength(c), pmt, z);  
                else
                    tempChunk = readSBX(sbxPath, sbxInfo, binLims(c,1), binLength(c), pmt, []); 
                end
                stackChan{pmt}(:,:,c) = mean(tempChunk, 3);
                waitbar((Nbin-c)/Nbin, w);
                %
            end
            delete(w);
            stackChan{pmt} = uint16(stackChan{pmt}); 
        end
        if any(edges > 0), stackChan{pmt} = stackChan{pmt}(edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:,:,:); end  % crop the data
        % resize  
        if scaleFactor ~= 1
            if Nscan > 1
                stackChan{pmt} = imresize3( stackChan{pmt}, [size(stackChan{pmt},1)/scaleFactor, size(stackChan{pmt},2)/scaleFactor, size(stackChan{pmt},3)] );
            else
                stackChan{pmt} = imresize( stackChan{pmt}, [size(stackChan{pmt},1)/scaleFactor, size(stackChan{pmt},2)/scaleFactor] );
            end
        end       
        % Save channel data to monochrome tif (optional)
        if monochrome
            if ~fileExists || overwrite
                if verbose, fprintf('\nWriting %s    ', pmtProjPath{pmt} ); end
                WriteTiff(stackChan{pmt}, pmtProjPath{pmt}); %pipe.io.writeTiff(stackChan{pmt}, chanProjPath{pmt});  %saveastiff( uint16(stackChan), chanProjPath{pmt}, RGBOpt ); %
            end
        end
    end
end 

% reverse color order to RGB
stackChan = stackChan([2,1]);
chanName = flip(pmtName);

%SWITCH ORDER TO RGB, MERGE CHANNELS AND (OPTIONAL) WRITE A SINGLE RGB TIFF
rgbPath = '';
if Npmt > 1
    % Write an RGB movie (optional)
    if RGBtoggle 
        rgbPath = sprintf('%s%s%s_Z%02d_RGB.tif', saveDir, saveName, sbxType, z);
        if (~exist(rgbPath,'file') || overwrite)
            if verbose, fprintf('\nWriting %s', rgbPath); end
            rgbTiff = cell(1,2);
            for pmt = 1:2
                if ~rescaleIntToggle
                    rgbTiff{pmt} = uint8(stackChan{pmt}/256);
                else
                    % Adjust data for uint8 RGB tiff if needed
                    chanLower = prctile(stackChan{pmt}(:), 1);
                    chanUpper = max(stackChan{pmt}(:)); %prctile(stackChan{chan}(:), 1);
                    fprintf('\nRescaling channel %s: [%i, %i] -> [0, 255]', chanName{pmt}, chanLower, chanUpper);
                    rgbTiff{pmt} = uint8(rescale(stackChan{pmt}, 0, 2^8-1, 'inputMin',chanLower)); % min(stackChan{chan}(:))
                end
            end
            rgbTiff = cat(4, rgbTiff{:});
            rgbTiff(end,end,end,3) = 0; % need to add a blank blue channel for tif format
            if verbose, fprintf('\nWriting %s    ', rgbPath ); end
            WriteTiff(rgbTiff, rgbPath); % pipe.io.writeTiff(rgbTiff, rgbPath);
        end
    end
    stackOut = cat(4, stackChan{:}); 
else
    chanInd = [2,1];
    stackOut = stackChan{chanInd(usePMT)}; 
end
if ~allowZeros,  stackOut(stackOut == 0) = NaN;  end % Set zeros to NaN to avoid averaging over blank frames elsewhere
if verbose, toc, end  
projPaths = [flip(pmtProjPath), {rgbPath}];

end