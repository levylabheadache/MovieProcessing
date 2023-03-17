function [edges, zUse] = SegmentCat3D(sbxInfo, varargin) % [zproj_interp, mask, trace_raw, trace_hp, segmentation_data, sbxPath] =
IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
%addParameter( IP, 'pmt', 1, @isnumeric )
addParameter( IP, 'chan', 'green', @ischar)
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'edges', [], @isnumeric )
addParameter( IP, 'zProj', [], @isnumeric )
addParameter( IP, 'projScans', [], @isnumeric )
addParameter( IP, 'censScans', [], @isnumeric )
addParameter( IP, 'blur', 2, @isnumeric ); %gaussian blur width for trace extraction
addParameter( IP, 'chunkSize', 930, @isnumeric )
addParameter( IP, 'minFoot', 50, @isnumeric )
addParameter( IP, 'hp_cutoff', 51, @isnumeric )
addParameter( IP, 'corrThresh', 75, @isnumeric )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxInfo, varargin{:} );
%pmt = IP.Results.pmt;
segName = IP.Results.name;
overwrite = IP.Results.overwrite;
% Determine which color to use
useChan = IP.Results.chan; %'green';
chanName = {'red','green'};
chanInd = find(strcmpi(useChan, chanName));
pmtName = flip(chanName);
pmt = find(strcmpi(useChan, pmtName));

% Determine file paths
if sbxInfo.Nplane == 1
    sbxSuff = 'sbxreg'; % 'sbx_affine'; %sprintf('%s%s.sbx_interp', sbxInfo.dir, sbxInfo.exptName );
else
    sbxSuff = 'sbx_interp';
end
paramsPath = sprintf('%s%s_seg_params%s.mat', sbxInfo.dir, sbxInfo.exptName, segName );
sbxPath = sprintf('%s%s.%s', sbxInfo.dir, sbxInfo.exptName, sbxSuff);
projPath = sprintf('%s%s_zproj%s.tif',sbxInfo.dir, sbxInfo.exptName, segName ); %sprintf('%s%s_zproj.mat',sbxInfo.dir, sbxInfo.exptName );
icaPath = sprintf('%s%s_2D_ica%s.mat', sbxInfo.dir, sbxInfo.exptName, segName); %strcat(fdir,filesep,fname,'_2D_ica.mat');
icaFiltPath = sprintf('%s%s_2D_ica_filt%s.mat', sbxInfo.dir, sbxInfo.exptName, segName);
prelim_path = sprintf('%s%s_seg_data_prelim%s.mat', sbxInfo.dir, sbxInfo.exptName, segName);
final_path = sprintf('%s%s_seg_data_final%s.mat', sbxInfo.dir, sbxInfo.exptName, segName);

% Save parameters used for segmentation
if ~exist(paramsPath,'file') || overwrite
    fprintf('\nSaving %s', paramsPath);
    save(paramsPath, 'IP');
else
    fprintf('\nLoading %s', paramsPath);
    load(paramsPath, 'IP');
end

zUse = IP.Results.zProj;
if isempty(zUse), zUse = 1:sbxInfo.Nplane; end

edges = IP.Results.edges;
if isempty(edges)
    warning('Edges are undefined, detecting with GetEdges');
    if sbxInfo.nchan == 1
        projStack = loadtiff( sprintf('%s%s_regProj.tif', sbxInfo.dir, sbxInfo.exptName) );
        edges = GetEdges(projStack(:,:,zUse(end)), 'minInt',1500, 'show',true); % true
    else
        projStack = loadtiff( sprintf('%s%s_regProj_%s.tif', sbxInfo.dir, sbxInfo.exptName, useChan) );
        edges = GetEdges(projStack(:,:,zUse(end)), 'minInt',1500, 'show',false); %
    end
end
projScans = IP.Results.projScans;
if isempty(projScans), projScans = 1:sbxInfo.totScan; end
censScans = IP.Results.censScans;
chunkSize = IP.Results.chunkSize;
hp_cutoff = IP.Results.hp_cutoff;
blurWidth = IP.Results.blur;
minFoot = IP.Results.minFoot;
corrThresh = IP.Results.corrThresh;
if sbxInfo.Nplane == 1, corrThresh = 0; end

% make sbx_interp if it doesn't exist
if sbxInfo.Nplane > 1 && (~exist(sbxPath, 'file')  || overwrite )
    InterpTemporalCat( sbxInfo, 'edges',edges, 'pmt',-1, 'chunkSize',6 );
end
%WriteSbxPlaneTif(sbxPath, sbxInfo, 8, 'chan','both', 'dir',strcat(sbxInfo.dir, 'Ztifs\Interp\'), 'name',sbxInfo.exptName, 'type','interp', 'verbose',true, 'overwrite',true); % , 'Nscan',64 , 'Nscan',100 , 'Nscan',4

% Run PCA/ICA on 2D projection
tic
if ~exist(icaFiltPath, 'file') || overwrite
    % Generate a projection, if one doesn't exist
    if sbxInfo.Nplane > 1
        fprintf('\nProjecting planes %i - %i', zUse(1), zUse(end) );
        zproj_interp = WriteSbxZproj(sbxPath, sbxInfo, 'z',zUse, 'chan','green', 'firstScan',projScans(1), 'Nscan',numel(projScans), 'type','zproj'); %
    else
        if sbxInfo.totScan < 20000
            binT = 1;
        else
            binT = 10;
        end
        zproj_interp = WriteSbxPlaneTif(sbxPath, sbxInfo, 1, 'chan','green', 'binT',binT, 'firstScan',projScans(1), 'Nscan',numel(projScans), 'type','zproj', 'overwrite',true ); % , 'name',sbxInfo.exptName
    end

    %{
    if ~exist(projPath, 'file') || overwrite
        if sbxInfo.Nplane > 1
            fprintf('\nProjecting planes %i - %i', zUse(1), zUse(end) );
            zproj_interp = WriteSbxZproj(sbxPath, sbxInfo, 'z',zUse, 'firstScan',projScans(1), 'Nscan',numel(projScans)-3, 'chan','green', 'type','interp', 'monochrome',false); %
            %zproj_interp = pipe.zproj(sbxPath, projScans(1), numel(projScans), pmt, zUse, 'mtype','.sbx_interp', 'registration',false);  %, 'write_tiff',true
            % Censor badly-registered scans
            if ~isempty(censScans)
                fprintf('\nCensoring scans %i - %i', censScans(1), censScans(end) );
                [~,projCensScans] = intersect(projScans, censScans ); % translate censored scans to the projScans frame
                zproj_interp(:,:,projCensScans) = [];
            end
            fprintf('\nWriting %s... ', projPath );
            WriteTiff(uint16(zproj_interp), projPath); %pipe.io.writeTiff(uint16(zproj_interp), projPath); toc
            toc
        elseif sbxInfo.totScan < 20000
            zproj_interp = WriteSbxPlaneTif(sbxPath, sbxInfo, 1, 'dir',sbxInfo.dir, 'name',sbxInfo.exptName, 'type','zproj', 'binT',1, 'verbose',true );
        else
            zproj_interp = WriteSbxPlaneTif(sbxPath, sbxInfo, 1, 'dir',sbxInfo.dir, 'name',sbxInfo.exptName, 'type','zproj', 'binT',10, 'verbose',true );
        end
        toc
    else
        fprintf('\nLoading %s...   ', projPath);
        tic
        zproj_interp = double(imread_big(projPath)); %double(pipe.io.readTiff(projPath)); % double();  loadtiff(projPath)
        toc
        if sbxInfo.Nplane > 1 % last scan of sbx_interp is bad
            tic
            zproj_interp(:,:,end) = [];
            toc
        end
    end
    %}

    % Run PCA/ICA on projection to generate initial 2D masks
    if ~exist(icaPath, 'file') || overwrite
        fprintf('\nPCA/ICA on projection');
        fix(clock)
        mask = PCAICA_2D_extract(zproj_interp, 'axons',true, 'hp_cutoff',hp_cutoff, 'edges',edges, 'save',icaPath);   % trace_raw, trace_hp
        fix(clock)
    else
        fprintf('\nLoading %s... ', icaPath);
        load(icaPath);
    end

    % Filter and break up masks
    ROI = FilterMasks(mask, 'minFoot',minFoot, 'edge',edges, 'break',false, 'show',true); % false
    ROIlabel = zeros( size(mask{1}));
    for m = 1:numel(ROI)
        ROIlabel( ROI(m).mask ) = m;
    end
    imshow( label2rgb(ROIlabel) );
    mask = {ROI.mask}; Nmask = numel(mask);

    % Extract mask traces from full data and high-pass filter them
    trace_raw = cell(1,Nmask);
    hp_filt = @(x)(x - movmedian(x, hp_cutoff));

    fprintf('\nExtracting traces');
    w = waitbar(0,'Extracting PCA mask traces...'); % w =
    tic
    for t = 1:size(zproj_interp,3)
        tempVol = zproj_interp(:,:,t);
        for m = 1:Nmask
            trace_raw{m}(1,t) = mean( tempVol(mask{m}(:)) );
        end
        waitbar(t/size(zproj_interp,3));
    end
    delete(w);
    %end
    trace_hp = cellfun(hp_filt, trace_raw, 'UniformOutput',false);
    fprintf('\nSaving extracted traces to %s... ', icaFiltPath);
    save(icaFiltPath, 'ROI', 'mask','Nmask', 'trace_raw', 'trace_hp', '-nocompression'); % , '-append'
else %if ~exist(prelim_path,'file') || overwrite
    fprintf('\nLoading %s... ', icaFiltPath);
    load(icaFiltPath);
    %FilterMasks(mask, 'minFoot',minFoot, 'edge',edges, 'break',false, 'show',true); % false
end
toc

tic
if ~exist(prelim_path,'file') || overwrite
    % Preliminary 3D segmentation
    fprintf('\nExtracting 3D masks (n = %i masks)', numel(mask));
    segmentation_data = Extract_3D_Mask(sbxPath, sbxInfo, mask, trace_hp, trace_raw, projScans, censScans, 'save',prelim_path, 'z',zUse, ...
        'blur',blurWidth, 'prctile_cutoff',corrThresh, 'chunkSize',chunkSize); %
    SegmentationROI3D( final_path, segmentation_data);  %
elseif ~exist(final_path, 'file') || overwrite
    % Curate
    fprintf('\nLoading %s... ', prelim_path);
    load( prelim_path ); toc % load segmentation_data variable
    SegmentationROI3D( final_path, segmentation_data);  %
end
toc
end


%{
    if sbxInfo.Nplane == 1 
        clearvars zproj_interp;
        if sbxInfo.totScan >= 20000
            Nchunk = 100;
        else
            Nchunk = 1;
        end
        % Break data into chunks and extract traces
        trace_raw_run = cell(Nmask, sbxInfo.Nrun);
        for r = 1:sbxInfo.Nrun
            % Read scans from sbx in chunks
            [chunkLims, Nchunk, chunkNscan] = MakeChunkLims(sbxInfo.scanLim(r), sbxInfo.scanLim(r+1)-1, sbxInfo.scanLim(r+1)-1, 'N',Nchunk ); % chunkSize, 
            trace_raw_chunk = cell(Nchunk, Nmask);  
            fprintf('\nLoading run %i, scans %i - %i, %i chunks\n', r, sbxInfo.scanLim(r), sbxInfo.scanLim(r+1)-1, Nchunk)
            tic 
            for c = 1:Nchunk
                fprintf('\nc = %i', c);
                zproj_chunk = double( readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkNscan(c), 1) ); 
                zproj_chunk = imgaussfilt(zproj_chunk, blurWidth); 
                zproj_chunk = reshape( zproj_chunk, [], chunkNscan(c) ); 
                for m = 1:Nmask % parfor is slower
                    trace_raw_chunk{c,m} = mask{m}(:)'*zproj_chunk;  % get raw mean trace from the movie chunk
                end
                toc
            end
            for m = 1:Nmask,  trace_raw_run{m,r} = cat(2, trace_raw_chunk{:,m});   end
        end
        clearvars zproj_chunk
        % Concatenate run traces and highpass filter
        for m = 1:Nmask
            trace_raw{m} = cat(2, trace_raw_run{m,:});
        end
        toc
%{
        clearvars zproj_interp;
        % Break data into chunks and extract traces
        Nchunk = 10;
        trace_raw_run = cell(Nmask, sbxInfo.Nrun);
        for r = 1:sbxInfo.Nrun
            [chunkLims, Nchunk, chunkNscan] = MakeChunkLims(sbxInfo.scanLim(r), sbxInfo.scanLim(r+1)-1, sbxInfo.scanLim(r+1)-1, 'N',Nchunk ); % chunkSize, 
            % Read scans from sbx in chunks
            fprintf('\nLoading run %i, scans %i - %i, %i chunks\n', r, sbxInfo.scanLim(r), sbxInfo.scanLim(r+1)-1, Nchunk)
            tic 
            zproj_chunk = cell(1,Nchunk); 
            for c = 1:Nchunk
                zproj_chunk{c} = double( readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkNscan(c), 1) ); %pipe.io.sbxRead(sbxPath, chunkLims(c,1), chunkLims(c,2)-chunkLims(c,1)+1, 1); 
                toc
            end
            %zproj_chunk = cellfun(@double, zproj_chunk, 'UniformOutput',false); toc
            zproj_chunk = cellfun(@imgaussfilt, zproj_chunk, repmat({blurWidth},1,Nchunk), 'UniformOutput',false); toc
            zproj_chunk = cellfun(@reshape, zproj_chunk, repmat({[]}, 1, Nchunk), num2cell(chunkNscan)', 'UniformOutput',false ); toc

            % get raw mean trace from the projection movie chunks
            fprintf('\nExtracting chunked raw traces\n');
            trace_raw_chunk = cell(Nchunk, Nmask); 
            tic
            for m = 1:Nmask % parfor is slower
               for c = 1:Nchunk % parfor is slower?
                   trace_raw_chunk{c,m} = mask{m}(:)'*zproj_chunk{c}; %reshape(zproj_chunk{c}, [], chunkNscan(c)
               end
               trace_raw_run{m,r} = cat(2, trace_raw_chunk{:,m});
            end
            toc
        end
        clearvars zproj_chunk
        % Concatenate run traces and highpass filter
        for m = 1:Nmask
            trace_raw{m} = cat(2, trace_raw_run{m,:});
        end
        toc
%}
    else
%}