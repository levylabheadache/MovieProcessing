function [regParam, projParam] = DefaultProcessingParams(dataSet)
% Registration parameters
regParam.refChan = '';
regParam.refScan = []; % which scans (in concatenated units) should be used to generate the reference image
regParam.refRun = 1; % which run should be used to generate the reference image (only used if refScan is empty)
regParam.chunkSize = 100; % register this many scans at once, to avoid loading too much data
regParam.edges = [80,80,40,40]; % cut this many pixels off from the [left, right, top, bottom] to avoid edges contaminating the registration
regParam.gamma = 1; % gamma correction parameter, used to suppress or accentuate brightest pixels
regParam.lowpass = 0; % use lowpass Gaussian spatial filtering to clean data before registration (pixels)
regParam.highpass = 0;% use highpass Gaussian spatial filtering to clean data before registration (pixels)
regParam.medFilter = [0,0,0]; % dimensions (pixels/scans) to use for median filtering to clean data during registration. set to zeros to skip med filtering
regParam.binXY = 2; % spatial downsampling factor - bigger number -> faster
regParam.avgT = 0; % 15;  temporal Gaussian filter width (not downsampling). scans, not seconds
regParam.avgTsigma = 0; % 5
regParam.binT = 1; % temporal downsampling - not currently implemented, keep at 1
regParam.histmatch = false; % true; % enable if registration struggles due to bleaching
regParam.method = ''; % use affine to get deformation, rigid just to align data
regParam.turboreg = true; %true
regParam.name = ''; % names can be used to distinguish different versions of registration

% Projection parameters
projParam.edge = [60,60,40,40];
projParam.type = 'max';
projParam.umPerPixel_target = 1; % target um/pix for projections. spatial downsampling 
switch lower(dataSet)
    case 'afferents'
        regParam.refChan = 'green'; 
        regParam.method = 'affine'; 
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'green'}; % 'red',
        projParam.vol = false;
    case 'affcsd'
        regParam.refChan = 'green';
        regParam.method = 'affine'; 
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'green'}; % 'red',
        projParam.vol = false;
    case 'astrocyte'
        regParam.refChan = 'red'; 
        regParam.method = 'translation'; 
        projParam.rate_target = 0.5; % Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    case 'macrophage'
        regParam.refChan = 'red';
        regParam.method = 'translation'; %translation
        projParam.rate_target = 1; % 1Hz
        projParam.color = {'red','green'};
        projParam.vol = false;
    case 'pollen'
        regParam.refChan = 'red';
        regParam.method = 'translation'; 
        projParam.rate_target = 1; % Hz
        projParam.color = {'red','green'}; % 'red',
        projParam.vol = true;
    case 'vasculature'
        regParam.refChan = 'red'; %red
        regParam.method = 'translation'; %'translation' 'affine'
        projParam.rate_target = 1; % Hz
        projParam.color = {'red'}; %'red','green'
        projParam.vol = false;
    otherwise
        regParam.refChan = 'green';
        regParam.method = 'translation';
        projParam.rate_target = 1; % Hz
        projParam.color = {'red','green'}; % 
        projParam.vol = false;
end
projParam.overwrite = false;
end