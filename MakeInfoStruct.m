function infoStruct = MakeInfoStruct(varargin) %  mainDir, mouse, exptDate, run, FOV
fSep = '\';
if nargin == 1
    if exist(varargin{1}, 'dir')
        runDir = varargin{1};
        dirSplit = split(runDir, '\');
        nameSplit = split( dirSplit{end-1}, '_');
        mouse = nameSplit{1};
        exptDate = nameSplit{2};
        runNum = str2double(nameSplit{3}(4:end));
        %strfind(nameSplit{3}, 'run');
    elseif exist(varargin{1}, 'file')
        [fDir, fName] = fileparts(varargin{1});
        runDir = strcat(fDir, fSep);
        nameSplit = split( fName, '_');
        mouse = nameSplit{1};
        exptDate = nameSplit{2};
        if numel(nameSplit{3}) > 3
            runNum = str2double(nameSplit{3}(4:end));
        else
            runNum = NaN;
        end
        infoStruct =  LoadSBXinfo(varargin{1}); %pipe.io.sbxInfo( varargin{1} );
    end
elseif nargin == 2
    runNum = varargin{2}; %expt.run;
    expt = varargin{1};
    if ~isstruct(expt), error('Single input must be the expt structure'), end
    mouse = expt.mouse;
    exptDate = expt.date;
    exptName = expt.name;
    exptDir = expt.dir;
    %fov = expt.fov;
    runDir = sprintf('%s%s_%s_run%i\\', exptDir, exptDate, mouse, runNum );
elseif nargin > 2
    mainDir = varargin{1};
    mouse = varargin{2};
    mouseDir = strcat(mainDir, mouse, fSep);
    exptDate = varargin{3};
    runNum = varargin{4};
    if isnumeric(exptDate), exptDate = num2str(exptDate); end
    if nargin > 4 && ~isempty(varargin{5})
        fov = num2str(varargin{5}); 
        exptDir = sprintf('%s%s_FOV%s\\', mouseDir, exptDate, fov); % sprintf('%s%s_FOV%i\\', mouseDir, exptDate, fov); 
        if ~exist(exptDir, 'dir')
            exptName = sprintf('%s_%s', mouse, exptDate);
            exptDir = sprintf('%s%s\\', mouseDir, exptDate); 
        else
            exptName = sprintf('%s_%s_FOV%s', mouse, exptDate, fov); 
        end
    else
        exptName = sprintf('%s_%s', mouse, exptDate); %sprintf('%s_%s', exptDate, mouse);
        exptDir = sprintf('%s%s\\', mouseDir, exptDate); % _%s , mouse
    end
    runDir = sprintf('%s%s_run%i\\', exptDir, exptName, runNum );
end
if ~exist(runDir, 'dir'), error('%s does not exist', runDir); end
%D:\2photon\CGRP02\201110_CGRP02\201110_CGRP02_run2
if ~exist('infoStruct', 'var')
    [~, sbxPath] = FileFinder(runDir, 'type','sbx' ); % 
    sbxPath = sbxPath{1};
    [fDir, fName] = fileparts( sbxPath );
    infoStruct = LoadSBXinfo(sbxPath); % pipe.io.sbxInfo( sbxPath ); 
    if isfield(infoStruct,'frame')
        infoStruct = rmfield(infoStruct, ["frame","event_id","line"]); % remove unused fields that sometimes cause problems
    end
end
sbxDir = strcat(fDir,fSep);
sbxPath = strcat(sbxDir, fName, '.sbx');
% Other useful info
infoStruct.mouse = mouse; 
infoStruct.exptDate = exptDate; 
infoStruct.run = runNum; 
infoStruct.Nrun = 1;
infoStruct.exptName = sprintf('%s_%s_run%i', infoStruct.mouse, infoStruct.exptDate, infoStruct.run );
infoStruct.fileName = fName;
infoStruct.dir = sbxDir;
infoStruct.path = sbxPath;

% Determine # of planes, frames, scans, framerate, digital zoom etc
if isfield( infoStruct, 'optotune_used' ) && ~infoStruct.optotune_used
    infoStruct.Nplane = 1;
else
    if isfield(infoStruct, 'otparam')
        infoStruct.Nplane = infoStruct.otparam(3);
    elseif isfield( infoStruct, 'otlevels' )
        infoStruct.Nplane = infoStruct.otlevels;
    elseif isfield( infoStruct, 'otwave' )
        if isempty(infoStruct.otwave)
            infoStruct.Nplane = 1;
        else
            infoStruct.Nplane = numel(infoStruct.otwave);
        end
    end
    if infoStruct.Nplane == 0, error('Nplane is 0!?'); end
end
if ~isfield( infoStruct, 'nframes' ) % || infoStruct.nframes = 
    infoStruct.nframes = infoStruct.config.frames;
end
infoStruct.totFrame = infoStruct.nframes;
infoStruct.Nscan = infoStruct.nframes/infoStruct.Nplane;
if rem(infoStruct.Nscan,1) ~= 0
    fprintf('\nNon-integer number of scans: rounding down to %i', floor(infoStruct.Nscan))
    infoStruct.Nscan = floor(infoStruct.Nscan); 
end
infoStruct.totScan = infoStruct.Nscan;
if ~isfield( infoStruct, 'framerate' )
    infoStruct.framerate = 15.49;
    fprintf('\nframerate field missing, assumed to be %2.2f', infoStruct.framerate );
end
if isfield(infoStruct.config, 'magnification_list')
    infoStruct.zoom = str2double(infoStruct.config.magnification_list(infoStruct.config.magnification,:));
    %fprintf('\nMagnification = %2.1f X', infoStruct.zoom);
else
    fprintf('\nMagnification not specified');
    infoStruct.zoom = NaN;
end
% Determine duration (seconds) and acquisition timestamp
infoStruct.duration = infoStruct.nframes/infoStruct.framerate;
metadata = dir( sbxPath );
endTime = datetime(metadata.date); %this is the time when the file was completed.
infoStruct.timestamp = endTime - seconds(infoStruct.duration); %this is the time when recording began

% Check if there is wheel data
%dataPath = FileFind(runDir, 'quadrature');

end

