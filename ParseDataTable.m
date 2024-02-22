function [expt, runInfo, varargout] = ParseDataTable(dataTable, x, dataCol, dataDir, varargin) % dataDir, projParam  dataTablePath, dataSet,
IP = inputParser;
addRequired( IP, 'dataTable', @iscell )
addRequired( IP, 'x', @isnumeric )
addRequired( IP, 'dataCol', @isstruct )
addRequired( IP, 'dataDir', @ischar )
addOptional( IP, 'regParams', struct(), @isstruct)
addOptional( IP, 'projParam', struct(), @isstruct)
parse( IP, dataTable, x, dataCol, dataDir, varargin{:} );  % , ROI
regParams = IP.Results.regParams;
projParam = IP.Results.projParam;

% Experiment metadata
expt.mouse = dataTable{x,dataCol.mouse};
expt.date = dataTable{x,dataCol.date}; 

if isnumeric(expt.date)
   expt.date = num2str(expt.date); 
end

if any(ismissing(dataTable{x,dataCol.FOV})) % check if an FOV # is specified
    expt.fov = 1;
    expt.dir = sprintf('%s%s\\%s\\', dataDir, expt.mouse, expt.date); % sprintf('%s%s\\%s_FOV%i_%s\\', dataDir, expt.mouse, expt.date, expt.fov, expt.mouse);
    expt.name = sprintf('%s_%s', expt.mouse, expt.date); 
else
    expt.fov = dataTable{x,dataCol.FOV};
    %if ischar(expt.fov), expt.fov = str2double(expt.fov); end
    if isnumeric(expt.fov)
       expt.fov = num2str(expt.fov); 
    end
    expt.dir = sprintf('%s%s\\%s_FOV%s\\', dataDir, expt.mouse, expt.date, expt.fov); % sprintf('%s%s\\%s_FOV%i_%s\\', dataDir, expt.mouse, expt.date, expt.fov, expt.mouse);
    expt.name = sprintf('%s_%s_FOV%s', expt.mouse, expt.date, expt.fov); 
end

if exist(expt.dir, 'dir')
    [~,exptPath] = FileFinder(expt.dir, 'contains','expt', 'type','mat');
    if ~isempty(exptPath)
        fprintf('\nLoading %s',exptPath{1})
        load(exptPath{1});
        expt = Expt;
    else
        fprintf('\n\nParsing %s\n', expt.name )
        runFolders = FileFinder(expt.dir, 'contains','run', 'type',0);
        [expt.runs, ~] = sort(cellfun(@GetRunNumber, runFolders, 'UniformOutput',true)', 'ascend'); %#ok<UDIM>  sortInd
        expt.Nruns = numel(expt.runs);

        % Check for Reference channel
        if ~isempty(dataCol.ref)
            expt.refChan = dataTable{x,dataCol.ref}; % determine which channel to use as reference for registration and concatenation
        else
            expt.refChan = 'green';
        end

        % Check for CSD runs
        missingData = cellfun(@all, cellfun(@ismissing, dataTable, 'UniformOutput',false), 'UniformOutput',true );
        if isempty(dataCol.csd) || missingData(x,dataCol.csd)
            expt.csd = NaN;  expt.Ncsd = 0;
        else
            expt.csd = dataTable{x,dataCol.csd};  
            expt.Ncsd = numel(expt.csd);
        end
        if isnan(expt.csd)
            expt.preRuns = expt.runs;
            expt.postRuns = [];
        else
            expt.preRuns = 1:expt.csd(1)-1;
            expt.postRuns = expt.csd(1):expt.Nruns; %expt.csd(1)+1 SCN 02/22/2024 - removed +1 to account CSD run as postRun (confirmed with Andy)
        end
        
        % Check for Vascular channel
        if ~isempty(dataCol.vascChan)
            expt.vascChan = dataTable{x,dataCol.vascChan}; 
        else
            expt.vascChan = 'green';
        end

    end

    % Get run metadata
    runInfo = GetRunInfo(expt, dataDir);

    % Concatenate runs metadata
    expt.Nrow = runInfo(1).sz(1); expt.Ncol = runInfo(1).sz(2);
    if all([runInfo.otlevels]==runInfo(1).otlevels) && all([runInfo.nchan]==runInfo(1).nchan)
        expt.Nplane = runInfo(1).Nplane; % otlevels
        expt.Nchan = runInfo(1).nchan;
    else
        error('All runs must have the same number of planes and channels')
    end
    expt.Nscan = floor([runInfo.nframes]/expt.Nplane); expt.totScan = sum(expt.Nscan); expt.totFrame = sum([runInfo.totFrame]);
    expt.scanLims = [0, cumsum(expt.Nscan)];
    expt.frameRate = runInfo(1).framerate;
    expt.scanRate = runInfo(1).framerate/expt.Nplane;
    expt.sbx.cat = strcat(expt.dir, expt.name, '.sbxcat '); %'.sbx_interp'
    expt.sbx.opt = strcat(expt.dir, expt.name, '.sbxopt ');
    expt.sbx.dft = strcat(expt.dir, expt.name, '.sbxdft ');
    expt.sbx.z = strcat(expt.dir, expt.name, '.sbxz '); %'.sbx_interp'
    expt.sbx.reg = strcat(expt.dir, expt.name, '.sbxreg ');
    expt.sbx.interp = strcat(expt.dir, expt.name, '.sbx_interp ');
    if isfield(runInfo(1).config, 'magnification_list')
        expt.zoom = str2double(runInfo(1).config.magnification_list(runInfo(1).config.magnification,:));
    elseif isfield(runInfo(1).config, 'magnification')
        warning('\nmagnification_list does not exist, using magnification field instead')
        %magList = [1, 1.2, 1.4, 1.7, 2, 2.4, 2.8, 3.4, 4, 4.8, 5.7, 6.7, 8];
        expt.zoom = runInfo(1).config.magnification; %magList()
    else
        expt.zoom = NaN;
    end
    expt.umPerPixel = (1/0.53)/expt.zoom;
    expt.refRun = floor(median(expt.runs));

    % Get regParam and projParam-related values from the spreadsheet, if they exist
    % Z projection plane sets
    if expt.Nplane > 1
        if isfield(dataCol, 'Zproj') && ~isempty(dataCol.Zproj) % expects the format:  'z1-z2, z3-z4, z5-z6'
            z_string = dataTable{x,dataCol.Zproj};
            comma_pos = [0, strfind(z_string, ','), length(z_string)+1];
            dash_pos = [strfind(z_string, '-')];
            Ndash = numel(dash_pos);
            z_proj = cell(1,Ndash);
            for d = 1:Ndash
                z_proj{d} = [str2double(z_string(comma_pos(d)+1:dash_pos(d)-1)), str2double(z_string(dash_pos(d)+1:comma_pos(d+1)-1))];
            end
            projParam.z = z_proj;
        end
    else
        projParam.z = {1};
    end

    % Check if a reference color for registration is specified
    if isfield(dataCol, 'ref') && ~isempty(dataCol.ref)
        regParams.refChan = dataTable{x,dataCol.ref};
    else
        regParams.refChan = 'green';
        warning('\nNo reference channel was specified, using %s by default', regParams.refChan);
    end
    
    % Crop this many pixels from left, right, top bottom for registration/projections
    if isfield(dataCol, 'edges') && ~isempty(dataCol.edges) && ~ismissing( dataTable{x,dataCol.edges})
        edge_str = dataTable{x,dataCol.edges}; % must be in the format of L, R, T, B
        comma_pos = [0, strfind(edge_str, ','), length(edge_str)+1]; % 
        for c = 1:length(comma_pos)-1
            regParams.edges(c) = str2double(edge_str(comma_pos(c)+1:comma_pos(c+1)-1)); 
        end
        projParam.edge = regParams.edges;
    end
    varargout{1} = regParams;
    varargout{2} = projParam;
else
    error('Experiment directory %s does not exist!', expt.dir)
end
end
%{
    if (expt.Nplane > 1) && ((isfield(dataCol, 'Ztop') && ~isempty(dataCol.Ztop)) || (isfield(dataCol, 'Zbot') && ~isempty(dataCol.Zbot)))
        zRange = sort([dataTable{x,dataCol.Ztop}, dataTable{x,dataCol.Zbot}]);
        projParam.z{1} = zRange(1):zRange(end);
    else
        projParam.z = {};
    end
%}