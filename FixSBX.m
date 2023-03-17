function sbxOutputInfo = FixSBX(sbxInputPath, sbxInputInfo, varargin)
checkInfo = @(x)(isstruct(x) || isempty(x));
IP = inputParser;
addRequired( IP, 'sbxInputPath', @ischar )
addRequired( IP, 'sbxInfo', checkInfo ) %  @isstruct
if isempty(sbxInputInfo)
    [pathDir, pathName, ~] = fileparts(sbxInputPath); % pathExt
    [~,infoPath] = FileFinder(pathDir, 'type','mat', 'criteria',@(x)(strcmp(x,pathName))); % 
    fprintf('\nNo info structure was provided. Loading %s', infoPath{1})
    sbxInputInfo = MakeInfoStruct( infoPath{1} );
end
addParameter( IP, 'z', 1:sbxInputInfo.Nplane, @isnumeric )
addParameter( IP, 'scans', 1:sbxInputInfo.Nscan, @isnumeric )
addParameter( IP, 'flip', false, @islogical )
addParameter( IP, 'proj', true, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxInputPath, sbxInputInfo, varargin{:} ); % mouse, exptDate,
zUse = IP.Results.z;
scanUse = IP.Results.scans;
flipZ = IP.Results.flip;
proj = IP.Results.proj;
overwrite = IP.Results.overwrite;

fprintf('\nFixing %s: flip = %i, scans %i - %i, z = %i - %i', sbxInputPath, flipZ, scanUse(1), scanUse(end), zUse(1), zUse(end))
sbxFixPath = strcat(sbxInputPath, 'fix'); %sbxInputPath; %
[~, chan] = DeterminePMT('both', sbxInputInfo); % usePMTind
if ~exist(sbxFixPath, 'file') || overwrite
    % Load the planar data, shift the planes
    tic;
    Nplane = numel(zUse);
    Nscan = numel(scanUse);
    planeData = cell(1,Nplane);
    k = 0;
    for z = zUse %1:sbxInfo.Nplane
        k = k+1;
        planeData{k} = WriteSbxPlaneTif(sbxInputPath, sbxInputInfo, z, 'verbose',true, 'chan',chan, 'firstScan',scanUse(1), 'Nscan',Nscan, 'monochrome',false, 'RGB',false, 'overwrite',true ); %, 'dir','', 'name',sbxInputInfo.fileName
    end
    toc

    % Adjust the metadata to reflect the new parameters (NOTE - other fields may need to be adjusted)
    sbxOutputInfo = sbxInputInfo;
    sbxOutputInfo.path = sbxFixPath;
    sbxOutputInfo.Nscan = Nscan;
    sbxOutputInfo.Nplane = Nplane;
    sbxOutputInfo.otlevels = Nplane;
    if flipZ,  sbxOutputInfo.otwave = flip(sbxOutputInfo.otwave);  end
    sbxOutputInfo.otwave = sbxInputInfo.otwave(zUse);
    sbxOutputInfo.nframes = sbxOutputInfo.Nscan*sbxOutputInfo.Nplane;
    sbxOutputInfo.max_idx = sbxOutputInfo.nframes-1;
    sbxOutputInfo.framerate = sbxInputInfo.framerate*(Nplane/sbxInputInfo.Nplane);

    % generate proj
    if proj
        projData = cellfun(@mean, planeData, repmat({3}, 1, numel(planeData)), 'UniformOutput',false);
        projData = cellfun(@squeeze, projData, 'UniformOutput',false);
        projData = cat(4, projData{:});
        projData = permute(projData, [1,2,4,3]); % [x,y,z,c]
        projData = circshift( projData, -1, 3);
        if flipZ, projData = flip(projData, 3); end
    end

    % reshape data to [chan x row x column x z x scan]
    rw = SbxWriter(sbxFixPath, sbxOutputInfo, '.sbxfix', true); % pipe.io.RegWriter(sbxFixPath, sbxInfo, '.sbxfix', true);
    if sbxOutputInfo.nchan > 1
        planeData = cat( 5, planeData{:} );
        planeData = flip( permute(planeData, [4,1,2,5,3]), 1); %  % switch from RGB to PMT chan order % permute(planeData, [4,1,2,5,3]);
        planeData = circshift( planeData, -1, 4);
        if flipZ, planeData = flip(planeData, 4); end
        planeData = reshape(planeData, [size(planeData,[1,2,3]), prod(size(planeData,[4,5]))] );
        tic
        fprintf('\nWriting sbxfix file');
        rw.write( planeData ); % rw.write(squeeze(uint16(tempScan)));
        toc
        rw.delete;
    else
        planeData = cat( 4, planeData{:} );
        planeData = circshift( planeData, -1, 4);
        planeData = flip(planeData, 4);
        planeData = permute(planeData, [1,2,4,3]); % permute(planeData, [5,1,2,4,3]); 
        planeData = reshape(planeData, [size(planeData,[1,2]), prod(size(planeData,[3,4]))] );
        fprintf('\nWriting %s', sbxFixPath);
        tic
        rw.write( planeData );
        toc
        rw.delete;
    end
    toc  
    fprintf('\nWrote %s', sbxFixPath)

    if Nplane ~= sbxInputInfo.Nplane || Nscan ~= sbxInputInfo.Nscan
        matFixPath = [sbxInputInfo.dir, sbxInputInfo.fileName, '.mat'];
        % Preserve the original metadata file
        if exist(matFixPath, 'file')
            origDir = [sbxInputInfo.dir,'Orig\']; mkdir(origDir)
            matOrigPath = [origDir, sbxInputInfo.fileName, '.mat']; %
            movefile(matFixPath, matOrigPath);
            fprintf('\nSaved %s', matOrigPath)
        end
        info = sbxOutputInfo;
        save(matFixPath, 'info');
        fprintf('\nSaved %s', matFixPath)
    end

    % Write mean projection tifs
    if proj
        % Set up proj file names
        projDir = sprintf('%s\\', fileparts(sbxInputPath));
        nameRoot = sprintf('%s_raw_meanProj_', sbxInputInfo.fileName);
        chanName = {'red','green'}; %pmtName = {'green','red'}; %
        % Save each channel as a 16bit, monochrome tif (optional), and rescale for RGB tif (optional)
        chanInd = [2,1];
        [usePMT, ~] = DeterminePMT('both', sbxInputInfo); % usePMT
        useChan = sort(chanInd(usePMT));
        Nchan = numel(useChan);
        rgbStack = zeros(sbxInputInfo.width, sbxInputInfo.height, sbxInputInfo.Nplane, 3);
        for chan = useChan
            chanProjPath = sprintf('%s%s%s.tif',  projDir, nameRoot, chanName{chan}); % sbxInputInfo.dir
            fprintf('\nWriting %s', chanProjPath);  % if verbose,  end
            if size(projData,4) > 1
                WriteTiff(uint16(projData(:,:,:,chan)), chanProjPath); % saveastiff(uint16(projData(:,:,:,chan)), chanProjPath{chan});
            else
                WriteTiff(uint16(projData), chanProjPath);
            end
            if Nchan > 1
                tempStack = projData(:,:,:,chan);
                chanLower = prctile(tempStack(:), 1);
                %chanUpper = max(tempStack(:)); %prctile(stackChan{chan}(:), 1);
                %if verbose, fprintf('\nRescaling %s channel: [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper); end
                rgbStack(:,:,:,chan) = rescale(tempStack, 0, 255, 'inputMin',chanLower); % min(stackChan{chan}(:))
            end
        end
        if Nchan > 1
            rgbPath = strcat(projDir, nameRoot,'RGB.tif');
            fprintf('\nWriting %s', rgbPath); %if verbose, end
            rgbStack = uint8(rgbStack);
            WriteTiff(rgbStack, rgbPath);
        end
    end
else
    fprintf('\n%s already exists', sbxFixPath);
end
end