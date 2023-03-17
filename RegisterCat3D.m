function sbxInfo = RegisterCat3D(sbxInfo, regParams, varargin) % mouse, exptDate
% Run scanbox data through stages of rigid 3D registration, z interpolation and affine registration. Optionally, write out tifs showing results from each stage.
IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'regParams', @isstruct )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'edges', [80,80,20,20], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'chunk', 20, @isnumeric ) %don't go over 20
%addParameter( IP, 'refChan', 'green', @ischar ) % for scanbox, 1 = green, 2 = red. -1 = both
%addParameter( IP, 'scale', 2, @isnumeric ) %
addParameter( IP, 'minInt', 1500, @isnumeric )
addParameter( IP, 'fix', false, @islogical )
addParameter( IP, 'flip', true, @islogical )
addParameter( IP, 'dewarp', 'affine', @ischar )
%addParameter( IP, 'binT', 1, @isnumeric )
parse( IP, sbxInfo, regParams, varargin{:} ); % mouse, exptDate,
%regParams = IP.Results.regParams;
overwrite = IP.Results.overwrite;
setEdges = IP.Results.edges;
fixSbx = IP.Results.fix;
flipZ = IP.Results.flip;
%binT = IP.Results.binT;
%scaleFactor = IP.Results.scale;
dewarpType = IP.Results.dewarp;
chunkSize = IP.Results.chunk;
minInt = IP.Results.minInt;
%regParams.refChan = IP.Results.regParams.refChan;  


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
if sbxInfo.Nplane > 1  && fixSbx
    if (~exist(sbxFixPath,'file') || overwrite)
        fprintf('\n   Correcting sbx z order... ');
        sbxInfo = FixSBX(sbxInputPath, sbxInfo, flipZ, overwrite);
    end
    sbxInputPath = sbxFixPath;
end

% Write tifs from the raw data
fprintf('\n   Writing raw projections... ');
WriteSbxProjection(sbxInputPath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','raw', 'overwrite',overwrite);  % writeChan rawProjPath,
if sbxInfo.Nplane == 1
    % Rigid, DFT-based alignment to most stable part of movie
    CorrectData2D(sbxInputPath, sbxInfo, regParams, shiftPath, sbxDftPath, 'overwrite',false, 'sigma',[3,3,5]);
    %WriteSbxPlaneTif(sbxInputPath, sbxInfo, 1, 'verbose',true, 'monochrome',true, 'RGB',false, 'type','raw', 'overwrite',overwrite, 'chan',regParams.refChan, 'edges',regParams.edges, 'binT',8); %
    WriteSbxProjection(sbxDftPath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','dft', 'overwrite',overwrite); % 
    %WriteSbxPlaneTif(sbxDftPath, sbxInfo, 1, 'chan',regParams.refChan, 'edges',regParams.edges, 'scale',regParams.binXY, 'monochrome',true, 'type','dft', 'overwrite',overwrite, 'verbose',true);
    sbxInfo.path = sbxDftPath;
else
    % lensing correction for neurolabware data
    if any(strcmpi(dewarpType, {'affine','rigid'}))
        if ~exist(sbxOptPath,'file') || overwrite
            fprintf('\n   Correcting Neurolabware data (%s)... ', dewarpType);
            GetOptotuneWarp(sbxInputPath, sbxInfo, 'chan',regParams.refChan, 'type',dewarpType, 'edges',setEdges, 'save',true);  % , 'show',true , 'scale',scaleFactor reg , 'firstRefScan',500
        end
        WriteSbxProjection(sbxOptPath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','opt', 'overwrite',overwrite);
        sbxInputPath = sbxOptPath;
    else
        fprintf('\ndewarp type not set to affine or rigid - dewarping skipped')
    end

    % 3D DFT shifts (rigid)
    [~,optProjPath] = FileFinder(fDir, 'contains',sprintf('opt_meanProj_%s', regParams.refChan));
    optProj = loadtiff( optProjPath{1} );
    optProjMean = mean(optProj, 3);
    if isempty(setEdges)
        optEdges = GetEdges( optProjMean, 'minInt',minInt, 'show',true ); % fprintf('edges = [%i, %i, %i, %i]\n', interpEdges )
    else
        optEdges = setEdges;
        ShowEdges(optEdges, optProjMean);
    end

    % Calculate rigid corrections
    tic
    if ~exist(shiftPath,'file') || overwrite
        fprintf('\n   Calculating  3D DFT shifts... ');
        CorrectData3D(sbxInputPath, sbxInfo, shiftPath, regParams.refChan, 'chunkSize',chunkSize, 'edges',optEdges, 'scale',regParams.binXY);
    end
    toc
    % make registered SBX file
    if ~exist(sbxDftPath,'file') || overwrite
        fprintf('\n   Writing sbxdft... ');
        MakeSbxDFT(sbxInputPath, sbxInfo, shiftPath, regParams.refChan, 'edges',optEdges, 'proj',true); % , 'zprojPath',zprojPath zproj_mean =
    end
    toc
    % write mean projection of rigid-corrected data
    WriteSbxProjection(sbxDftPath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','dft', 'overwrite',overwrite);

    % z interpolation
    [~,dftProjPath] = FileFinder(fDir, 'contains',sprintf('dft_meanProj_%s', regParams.refChan));
    dftProj = loadtiff( dftProjPath{1} );
    dftProjMean = mean(dftProj, 3);
    if isempty(setEdges)
        dftEdges = GetEdges( dftProjMean, 'minInt',minInt, 'show',true ); % fprintf('edges = [%i, %i, %i, %i]\n', interpEdges )
    else
        dftEdges = setEdges;
        ShowEdges(dftEdges, dftProjMean);
    end

    if ~exist(sbxZpath,'file') || overwrite
        if ~exist(interpPath,'file')
            InterpZ(sbxDftPath, sbxInfo, interpPath, 'scale',regParams.binXY, 'edges',dftEdges, 'chunkSize',chunkSize); %DFT_reg_z_interp(sbxDftPath, interpMatPath, refChan, regParams.binXY, Nchunk, 'optotune',false, 'edges',dftEdges);
        end
        MakeSbxZ(sbxDftPath, sbxInfo, interpPath); % SBX_z_interp(sbxDftPath, interpMatPath);
    end
    WriteSbxProjection(sbxZpath, sbxInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','Z', 'overwrite',overwrite);
end
end