function [info, expt] = ConcatenateRunInfo(expt, runInfo, varargin)

IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'runInfo', @isstruct )
addParameter( IP, 'suffix', 'sbx', @ischar )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, runInfo, varargin{:} ); 
suffix = IP.Results.suffix;
overwrite = IP.Results.overwrite;

%Generate a metadata structure for concatenated runs
catInfoPath = sprintf('%s%s.mat', expt.dir, expt.name);
if ~exist(catInfoPath, 'file') || overwrite
    % Check that key parameters are constant between runs, otherwise produce an error
    if numel(unique([runInfo.scanmode])) > 1, error('Scanmode is inconsistent between runs'); end
    if numel(unique([runInfo.Nplane])) > 1, error('Nplane is inconsistent between runs'); end
    if numel(unique([runInfo.nchan])) > 1, error('Nchan is inconsistent between runs'); end
    
    info = runInfo(1);
    if ~isfield(info, 'calibration'), info.calibration = []; end
    if ~isfield(info, 'objective'), info.objective = ''; end
    info.fid = []; %fopen( )
    info.run = [runInfo.run];
    info.Nrun = numel(info.run);
    info.runDir = {runInfo.dir}';
    info.dir = expt.dir; 
    info.exptName = expt.name; %exptName; 
    info.fileName = info.exptName;
    %info.path = sprintf('%s%s.%s', info.dir, info.exptName, suffix); %expt.sbx; %
    info.path = sprintf('%s%s.%s', info.dir, info.exptName, suffix); %expt.sbx; %  expt.fov _FOV%i
    info.runFrames = [runInfo.nframes]; 
    %info.Nplane = runInfo(1).Nplane;
    info.nframes = sum(info.runFrames); % various pipe functions expect this field
    info.max_idx = info.nframes - 1;
    info.Nscan = [runInfo.Nscan]; 
    info.totScan = sum(info.Nscan);
    info.totFrame = sum(info.runFrames);
    info.scanLim = [1, cumsum(info.Nscan)+1];
    info.duration = sum([runInfo.duration]);
    info.timestamp = [runInfo.timestamp];

    fprintf('\nSaving %s', catInfoPath );
    save( catInfoPath, 'info' );
else
    fprintf('\nLoading %s', catInfoPath );
    load(catInfoPath, 'info')
end

expt.Nx = info.sz(1); 
expt.Ny = info.sz(2); 
expt.Nplane = info.Nplane; 
expt.Nchan = info.nchan;
expt.Nscan = info.Nscan; %floor([runInfo{x}.nframes]/expt.Nplane); 
expt.totScan = sum(expt.Nscan); 
expt.scanLims = [0, cumsum(expt.Nscan)];
expt.scanRate = info.framerate/expt.Nplane;
if isfield(info.config, 'magnification_list')
    expt.zoom = str2double(info.config.magnification_list(info.config.magnification,:));
else
    fprintf('\nMagnification unclear, assumed to be 2.4 X');
    expt.zoom = 2.4;
end
expt.umPerPixel = (1/0.53)/expt.zoom;

end