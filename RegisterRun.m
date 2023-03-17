function sbxInfo = RegisterRun(sbxInfo, regParams, varargin) % mouse, exptDate
% Run scanbox data through stages of optotune correction, rigid registration, and z interpolation. Also, write out tifs summarizing results from each stage.
IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'regParams', @isstruct )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'chunk', 20, @isnumeric ) %don't go over 20
addParameter( IP, 'fix', false, @islogical )
addParameter( IP, 'flip', true, @islogical )
addParameter( IP, 'dewarp', 'affine', @ischar )
addParameter( IP, 'interp', true, @islogical )
parse( IP, sbxInfo, regParams, varargin{:} ); % mouse, exptDate,
overwrite = IP.Results.overwrite;
setEdges = regParams.edges; %IP.Results.edges;
fixSbx = IP.Results.fix;
flipZ = IP.Results.flip;
dewarpType = IP.Results.dewarp;
chunkSize = IP.Results.chunk;
toggleInterp = IP.Results.interp;
refChanInd = find(strcmpi({'red','green'}, regParams.refChan));

% Determine file paths and names
sbxInputPath = sbxInfo.path;
[fDir, fName] = fileparts(sbxInputPath);
pathTemplate = strcat( fDir, '\', fName );
sbxFixPath = [pathTemplate, '.sbxfix'];
sbxOptPath = [pathTemplate, '.sbxopt'];
sbxDftPath = [pathTemplate, '.sbxdft'];
sbxZpath = [pathTemplate, '.sbxz'];
shiftPath = strcat(pathTemplate, '_dftshifts.mat');
interpPath = strcat(sbxInfo.dir,sbxInfo.exptName,'_zinterp.mat');

if rem(sbxInfo.totScan, chunkSize) ~= 0
    chunkSize = max(factor(sbxInfo.totScan));
    fprintf('Non-integer number of chunks, changed chunk size to %i', chunkSize);
end
fprintf('\nProcessing %s: %i frames, %i planes, %i scans. %i channel\n', sbxInfo.exptName, sbxInfo.nframes, sbxInfo.Nplane, sbxInfo.totScan, sbxInfo.nchan)
if sbxInfo.Nplane > 1 && fixSbx
    if (~exist(sbxFixPath,'file') || overwrite)
        fprintf('\n   Correcting sbx z order... ');
        sbxInfo = FixSBX(sbxInputPath, sbxInfo, flipZ, overwrite);
    else
        sbxInfo.path = sbxFixPath;
    end
end
if exist(sbxFixPath, 'file')
    %sbxInputPath = sbxFixPath;
    sbxInfo.path = sbxFixPath;
end


% Write tifs from the raw data
rawProjMean = WriteSbxProjection(sbxInfo.path, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','raw', 'overwrite',overwrite);  % writeChan rawProjPath,
if isempty(regParams.edges)
    %fprintf('\n   Writing raw projections... ');
    regParams.edges = GetEdges3D(rawProjMean(:,:,:,refChanInd), 'rate_max',1, 'allow_zero',false, 'show',true);
end

if sbxInfo.Nplane == 1
    % Rigid, DFT-based alignment to most stable part of movie
    % If the reference scans have been provided, use them to generate a reference image
    refImage = WriteSbxPlaneTif(sbxPath, sbxInfo, 1, 'chan',regParams.refChan, 'edges',regParams.edges, 'monochrome',false, 'RGB',false, 'overwrite',true, 'verbose',verbose);


    CorrectData2D(sbxInfo.path, sbxInfo, regParams, shiftPath, sbxDftPath, 'overwrite',false, 'sigma',[3,3,5]);
    sbxInfo.path = sbxDftPath;
    %WriteSbxPlaneTif(sbxInputPath, sbxInfo, 1, 'verbose',true, 'monochrome',true, 'RGB',false, 'type','raw', 'overwrite',overwrite, 'chan',regParams.refChan, 'edges',regParams.edges, 'binT',8); %
    WriteSbxProjection(sbxInfo.path, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','dft', 'overwrite',overwrite); %
    %WriteSbxPlaneTif(sbxDftPath, sbxInfo, 1, 'chan',regParams.refChan, 'edges',regParams.edges, 'scale',regParams.binXY, 'monochrome',true, 'type','dft', 'overwrite',overwrite, 'verbose',true);
else
    if (~exist(sbxZpath, 'file') || ~exist(interpPath, 'file')) || overwrite
        if (~exist(sbxDftPath, 'file') || ~exist(shiftPath, 'file')) || overwrite
            % lensing correction for neurolabware data
            if any(strcmpi(dewarpType, {'affine','rigid'}))
                if ~exist(sbxOptPath,'file') || overwrite
                    fprintf('\n   Correcting Neurolabware data (%s)... ', dewarpType);
                    GetOptotuneWarp(sbxInfo.path, sbxInfo, 'firstRefScan',regParams.refScan(1), 'Nref',numel(regParams.refScan), 'chan',regParams.refChan, 'type',dewarpType, 'edges',setEdges, 'save',true);  % , 'show',true , 'scale',scaleFactor reg , 'firstRefScan',500
                end
                sbxInfo.path = sbxOptPath;
            else
                fprintf('\ndewarp type not set to affine or rigid - dewarping skipped')
            end
            % 3D DFT shifts (rigid)
            optProjMean = WriteSbxProjection(sbxOptPath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','opt', 'overwrite',overwrite);
            if isempty(setEdges)
                %optEdges = GetEdges( optProjMean, 'minInt',minInt, 'show',true ); % fprintf('edges = [%i, %i, %i, %i]\n', interpEdges )
                optEdges = GetEdges3D(optProjMean(:,:,:,refChanInd), 'rate_max',1, 'allow_zero',false, 'show',true); 
            else
                optEdges = setEdges;
                ShowEdges(optEdges, optProjMean(:,:,end));
            end
        else
            sbxInfo.path = sbxOptPath;
        end

        % Calculate rigid corrections
        if ~exist(shiftPath,'file') || overwrite
            fprintf('\n   Calculating  3D DFT shifts... ');  tic
            CorrectData3D(sbxInfo.path, sbxInfo, shiftPath, regParams.refChan, 'chunkSize',chunkSize, 'edges',optEdges, 'scale',regParams.binXY);
            toc
        end
        
        % make sbxdft file
        if ~exist(sbxDftPath,'file') || overwrite
            fprintf('\n   Writing sbxdft... '); tic
            MakeSbxDFT(sbxInfo.path, sbxInfo, shiftPath, regParams.refChan, 'edges',optEdges, 'proj',true); % , 'zprojPath',zprojPath zproj_mean =
            toc
        end
        sbxInfo.path = sbxDftPath;
        % write mean projection of rigid-corrected data
        dftProjMean = WriteSbxProjection(sbxDftPath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','dft', 'overwrite',overwrite);
        
        % z interpolation (optional)
        if toggleInterp
            if isempty(setEdges)
                %dftEdges = GetEdges( dftProjMean, 'minInt',minInt, 'show',true ); % fprintf('edges = [%i, %i, %i, %i]\n', interpEdges )
                dftEdges = GetEdges3D(dftProjMean(:,:,:,refChanInd), 'rate_max',1, 'allow_zero',false, 'show',true);
            else
                dftEdges = setEdges;
                ShowEdges(dftEdges, dftProjMean(:,:,end));
            end

            if ~exist(sbxZpath,'file') || overwrite
                if ~exist(interpPath,'file')
                    InterpZ(sbxDftPath, sbxInfo, interpPath, 'scale',regParams.binXY, 'edges',dftEdges, 'chunkSize',chunkSize); %DFT_reg_z_interp(sbxDftPath, interpMatPath, refChan, regParams.binXY, Nchunk, 'optotune',false, 'edges',dftEdges);
                end
                MakeSbxZ(sbxDftPath, sbxInfo, interpPath); % SBX_z_interp(sbxDftPath, interpMatPath);
            end
            sbxInfo.path = sbxZpath;
            WriteSbxProjection(sbxInfo.path, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','Z', 'overwrite',overwrite);
        end
    else
        sbxInfo.path = sbxZpath;
    end
end
close all;
end