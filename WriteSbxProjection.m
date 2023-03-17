function [projStack, rgbStack, projPath] = WriteSbxProjection(sbxPath, sbxInfo, varargin) % projPath, 
% Extract data from an SBX file, and then perform a projection (mean or max) across frames
checkInfo = @(x)(isstruct(x) || isempty(x));
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', checkInfo )
%addRequired( IP, 'sbxInfo', checkInfo ) % @isstruct
if isempty(sbxInfo)
    [pathDir, pathName, ~] = fileparts(sbxPath); % pathExt
    [~,infoPath] = FileFinder(pathDir, 'type','mat', 'criteria',@(x)(strcmp(x,pathName))); % 
    fprintf('\nNo info structure was provided. Loading %s', infoPath{1})
    sbxInfo = MakeInfoStruct( infoPath{1} );
end
%addOptional( IP, 'projPath', '', @ischar)
addParameter( IP, 'chan', 'both', @ischar ) % 'green', 'red' or 'both'. for scanbox, PMT 1 = green, PMT2 = red. -1 = both. RGB for channels 
addParameter( IP, 'z', [], @isnumeric ) 
addParameter( IP, 'projType', 'mean', @ischar)
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'Nscan', -1, @isnumeric )
addParameter( IP, 'scanMax', 10000, @isnumeric )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric ) % [left, right, top, bottom]
addParameter( IP, 'scale', 1, @isnumeric ) %
addParameter( IP, 'binT', 1, @isnumeric ) 
addParameter( IP, 'dir', '', @ischar ) 
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'type', '', @ischar ) 
addParameter( IP, 'monochrome', false, @islogical ) % write a monochrome tif for each plane?
addParameter( IP, 'RGB', false, @islogical ) % write an RGB tif for each plane?
addParameter( IP, 'rescale', true, @islogical )
addParameter( IP, 'verbose', true, @islogical );
addParameter( IP, 'overwrite', false, @islogical );
parse( IP, sbxPath, sbxInfo, varargin{:} );  % projPath, 
% Projection parameters
projChan = IP.Results.chan;  
projType = IP.Results.projType;
firstScan = IP.Results.firstScan;
edges = IP.Results.edges;
scaleFactor = IP.Results.scale;
binT = IP.Results.binT;
Nscan = IP.Results.Nscan;
scanMax = IP.Results.scanMax;
if Nscan == -1, Nscan = sbxInfo.totScan-firstScan+1; end
if Nscan > scanMax && binT == 1
    binT = round(Nscan/scanMax);  
    fprintf('\nNscan exceeds %i: setting binning to %i', scanMax, binT);
end  
zSet = IP.Results.z;
if isempty(zSet), zSet = 1:sbxInfo.Nplane; end
chanName = {'red','green'}; %pmtName = {'green','red'}; %
chanInd = [2,1];
[usePMT, usePMTname] = DeterminePMT(projChan, sbxInfo); % usePMT
useChan = sort(chanInd(usePMT));
Nchan = numel(useChan);
% Determine output file path and name
saveDir = IP.Results.dir; 
if ~isempty(saveDir), mkdir(saveDir); end
saveName = IP.Results.name;
if isempty(saveName), [~,saveName] = fileparts(sbxPath); end
if isempty(saveDir), saveDir = strcat(fileparts(sbxPath), '\'); end
sbxType = IP.Results.type;
if ~isempty(sbxType)
    nameRoot = sprintf('%s_%s_%sProj_', saveName, sbxType, projType); 
else
    nameRoot = sprintf('%s_%sProj_', saveName, projType);
end
% Determine the outputs
overwrite = IP.Results.overwrite;
monochrome = IP.Results.monochrome;
RGB = IP.Results.RGB;
if Nchan < 2, RGB = false; end
rescaleIntToggle = IP.Results.rescale;
verbose = IP.Results.verbose;
if scaleFactor == 1
    projStack = zeros(sbxInfo.width, sbxInfo.height, sbxInfo.Nplane, 2); % color order: red, green
    rgbStack = zeros(sbxInfo.width, sbxInfo.height, sbxInfo.Nplane, 3);
end
%projPath = cell(1,3);
if nargout == 0 && ~(monochrome || RGB) 
    if verbose, fprintf('\nNo output was requested - skipping\n'); end
else
    % Determine which channel(s) to use, create paths to monochrome projections and check if they already exist
    chanProjPath = cell(1,2); chanProjExists = false(1,2);
    for chan = useChan 
        chanProjPath{chan} = sprintf('%s%s%s.tif', saveDir, nameRoot, chanName{chan}); % chanName
        if exist(chanProjPath{chan},'file') && ~overwrite
            chanProjExists(chan) = true;
            if verbose, fprintf('\nLoading %s', chanProjPath{chan}); end
            projStack(:,:,:,chan) = loadtiff(chanProjPath{chan}); % projStack
        end
    end
    % If necessary, calculate the projection
    if ~all(chanProjExists(useChan)) || overwrite
        maxZ = @(x)(max(x,[],3));
        % Get each plane, crop, resize, and mean project (can't necessarily get the full data at once due to memory constraints)
        if verbose, fprintf('\nCalculating %s projection', projType); end
        for Z = flip(1:numel(zSet)) % flip(zSet) %flip(1:sbxInfo.otlevels)
            if verbose, fprintf('\nZ = %i', Z); end
            [~, tempChan] = WriteSbxPlaneTif(sbxPath, sbxInfo, zSet(Z), 'verbose',verbose, 'dir',saveDir, 'name',saveName, 'overwrite',overwrite, ...
                'edges',edges, 'scale',scaleFactor, 'firstScan',firstScan, 'Nscan',Nscan, 'chan',usePMTname, 'zeros',true, 'binT',binT, 'RGB',false, 'monochrome',false ); %
            if strcmpi(projType,'mean')
                tempProj = cellfun(@mean, tempChan, {3,3}, 'UniformOutput',false);
            else
                tempProj = cellfun(maxZ, tempChan, 'UniformOutput',false);
            end
            for chan = useChan
                projStack(:,:,Z,chan) = tempProj{chan};
            end
        end
        projStack = projStack(:,:,:,useChan); %usePMT
        % Save the results to monochrome and/or RGB tifs, as requested
        if Nchan == 1
            if monochrome
                if verbose, fprintf('\nWriting %s', chanProjPath{useChan}); end
                WriteTiff(uint16(projStack(:,:,:,1)), chanProjPath{useChan});  %saveastiff(uint16(projStack(:,:,:,1)), chanProjPath{useChan});
            end
        else
            % Save each channel as a 16bit, monochrome tif (optional), and rescale for RGB tif (optional)
            for chan = useChan
                if monochrome && ~all(chanProjExists(useChan))
                    if verbose, fprintf('\nWriting %s', chanProjPath{chan}); end
                    WriteTiff(uint16(projStack(:,:,:,chan)), chanProjPath{chan}); % saveastiff(uint16(projStack(:,:,:,chan)), chanProjPath{chan});
                end
                if RGB
                    if ~rescaleIntToggle
                        rgbStack(:,:,:,chan) = im2uint8(projStack(:,:,:,chan)); % /255
                    else
                        tempStack = projStack(:,:,:,chan);
                        chanLower = prctile(tempStack(:), 1);
                        %chanUpper = max(tempStack(:)); %prctile(stackChan{chan}(:), 1);
                        %if verbose, fprintf('\nRescaling %s channel: [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper); end
                        rgbStack(:,:,:,chan) = rescale(tempStack, 0, 255, 'inputMin',chanLower); % min(stackChan{chan}(:))
                    end
                end
            end
        end
    else
        projStack = projStack(:,:,:,useChan); %usePMT
    end
    % Save a combined RGB tif (optional)
    rgbPath = sprintf('%s%sRGB.tif', saveDir, nameRoot);
    if RGB && (~exist(rgbPath, 'file') || overwrite || ~all(chanProjExists(useChan)))
        if verbose, fprintf('\nWriting %s', rgbPath); end
        rgbStack = uint8(rgbStack);
        WriteTiff(rgbStack, rgbPath);
    else
        rgbPath = '';
    end
end
projPath = [chanProjPath, {rgbPath}];
end