function [reg_transforms, varargout] = RegisterSBX(mov_path, sbxInfo, refImage, params, varargin)
%   NOTE: hardcoded path to ImageJ.
IP = inputParser;
addRequired( IP, 'mov_path', @ischar )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'refImage', @isnumeric )
addRequired( IP, 'params', @isstruct )
addOptional(IP, 'z', 1, @isnumeric)  % Which optotune level to align
addParameter(IP, 'method', 'affine', @ischar )
addParameter(IP, 'firstScan', 1)  % The frame to start reading from
addParameter(IP, 'Nscan', 1)  % Number of frames to read 
addParameter(IP, 'overwrite', false, @islogical)
addParameter(IP, 'intTif', false, @islogical)
addParameter(IP, 'finalTif', false, @islogical)
addParameter(IP, 'verbose', false, @islogical)
addParameter(IP, 'dir', '', @ischar)
addParameter(IP, 'name', '', @ischar)
parse( IP, mov_path, sbxInfo, refImage, params, varargin{:} ); 
z = IP.Results.z;
firstScan = IP.Results.firstScan;
Nscan = IP.Results.Nscan;
lastScan = firstScan + Nscan - 1;
saveName = IP.Results.name;
if isempty(params.name)
    nameStr = '';
else
    nameStr = ['_',params.name];
    saveName = [saveName, nameStr];
end
writeIntermediateTif = IP.Results.intTif;
writeFinalTif = IP.Results.finalTif;
verbose = IP.Results.verbose;
overwrite = IP.Results.overwrite;
RegTempDir = sprintf('%sRegTemp\\', sbxInfo.dir); [~,~] = mkdir(RegTempDir);
chunkTempDir = sprintf('%sChunks\\', RegTempDir); [~,~] = mkdir(chunkTempDir);
chanName = {'green','red'};
params.refChanInd = find( strcmpi(params.refChan, chanName));

%tic
% Crop the reference image
rawRef = uint16(refImage(params.edges(3)+1:end-params.edges(4), params.edges(1)+1:end-params.edges(2)));
preRef = rawRef;

% Break the data into chunk of scans to avoid memory issues and improve TurboReg performance
if params.chunkSize == 0, params.chunkSize = Nscan; end % +1   sbxInfo.totScan
if params.chunkSize > Nscan, params.chunkSize = Nscan; end
[chunkScan, Nchunk, chunkLength] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',params.chunkSize, 'allowPartial',true );

if nargout > 1 || writeFinalTif
    unregData = zeros(size(rawRef,1), size(rawRef,2), numel(firstScan:lastScan) );
    regData = zeros(size(rawRef,1), size(rawRef,2), numel(firstScan:lastScan) );
end

% Set up registration method
if strcmpi(params.method, 'affine') && params.turboreg
    % Run imagej and prepare turboreg
    if verbose, fprintf('\nRunning ImageJ  '); end
    RunIJ();  % Hardcoded path to ImageJ
    trhl = TurboRegHL_();
else
    %[~, metric] = imregconfig('monomodal'); % optimizer
    %optimizer = registration.optimizer.OnePlusOneEvolutionary;
    optimizer = registration.optimizer.RegularStepGradientDescent;
    optimizer.MaximumIterations = 1000; % increase to possibly improve performance, at the expense of time (default = 100)
    optimizer.MinimumStepLength = 1.0000e-06; % decrease to possibly improve performance, at the expense of time (default = 1.000000e-05)
    optimizer.MaximumStepLength = 0.001; % decrease to possibly improve performance, at the expense of time (default = 0.0625)
    optimizer.RelaxationFactor = 0.75; % increase to improve performance on noisy data, at the expense of time (default = 0.5)
    metric = registration.metric.MeanSquares; % registration.metric.MattesMutualInformation is slower and doesn't seem to perform better
    targets = '';
    targetstr = '';
    szstr = '';
end
alignstr = '';

% Set up temporal averaging
weightVector = ones(1, params.avgT);
if params.avgT > 0
    if ( rem(params.avgT,2) == 0 )  
        params.avgT = params.avgT+1;
        warning('params.avgT should be set to an odd number - setting to %i', params.avgT);  
    end 
    avgTpad = round(params.avgT/2) - 1; % how many scans before/after the scan do we average over?
    if params.avgTsigma > 0
        tGauss = -avgTpad:avgTpad; 
        weightVector = exp(-(tGauss).^2/(2*params.avgTsigma^2));  
        weightVector = double(weightVector/sum(weightVector)); %plot(tGauss, weightVector)       
    end
else
    avgTpad = 0;
end

% Run the registration procedure, chunk by chunk
disp(params);
%{
if Nscan > 1
    w = waitbar(0,sprintf('Performing registration, %s', saveName));
    w.Children.Title.Interpreter = 'none';
end
%}
reg_transforms = repmat( {affine2d}, 1, Nscan ); % reg_transforms = cell(1, sbxInfo.totScan);
K = 0; % counter to keep track of regData frames
for c = 1:Nchunk
    chunkSavePath = sprintf('%sz%02i%s_chunk%04i.mat', chunkTempDir, z, nameStr, c);
    chunkScans = chunkScan(c,1):chunkScan(c,2);
    NchunkScans = chunkLength(c); %size(chunkScans, 3); 
    padScans = chunkScan(c,1)-avgTpad:chunkScan(c,2)+avgTpad; %
    padScans(padScans < 1) = []; padScans(padScans > sbxInfo.totScan) = [];
    NpadScans = numel(padScans);
    if ~exist(chunkSavePath, 'file') || overwrite
        % LOAD DATA
        if verbose, fprintf('\nchunk %i / %i:  Loading %s (z = %i, scans %i - %i)', c, Nchunk, mov_path, z, chunkScan(c,1), chunkScan(c,2)  ); end
        if sbxInfo.Nplane > 1
            rawData = readSBX(mov_path, sbxInfo, padScans(1), numel(padScans), params.refChanInd, z); 
        else
            rawData = readSBX(mov_path, sbxInfo, padScans(1), numel(padScans), params.refChanInd, []);
        end
        if numel(rawData) == 0; error('Blank rawData!'); end
 
        % PRE-PROCESS DATA AND REFERENCE IMAGE (note, reference image gets unnecessarily processed for each chunk, consider preprocessing once outside chunk loop)
        % Crop edges
        %if verbose, fprintf('\nCropping edges: [L, R, T, B] = [%i, %i, %i, %i]', params.edges(1), params.edges(2), params.edges(3), params.edges(4) ); end
        rawData = rawData(params.edges(3)+1:end-params.edges(4), params.edges(1)+1:end-params.edges(2), :);
        preData = rawData;
        %figure; subplot(1,2,1); imshow(rawRef,[]); subplot(1,2,2); imshow(preRef,[]); 

        % Gamma correction (suppress or accentuate very bright points)
        if params.gamma ~= 1
            if c == 1
                % Gamma correct the reference image
                gammaRef = imadjust(preRef,[],[],params.gamma);
                preRef = gammaRef;
                %{
                figure; 
                subplot(2,2,1); imshow(preRef, []); title('Reference')
                subplot(2,2,2); imshow(gammaRef, []); title(sprintf('Gamma = %2.2f', params.gamma))
                impixelinfo;
                subplot(2,2,3); imhist(preRef); subplot(2,2,4); imhist(gammaRef);
                %}
            end
            % Gamma correct the current data chunk
            for fr = 1:size(preData,3)
                preData(:,:,fr) = imadjust(preData(:,:,fr),[],[],params.gamma);
            end
        end
        
        % Low-pass filter 
        if params.lowpass > 0
            if verbose, fprintf('\nLowpass = %2.1f pix   ', params.lowpass ); end
            if c == 1
                preRef = imgaussfilt(preRef, params.lowpass);
            end
            preData = imgaussfilt(preData, params.lowpass);
        end

        % High-pass filter
        if params.highpass > 0
            if verbose, fprintf('\nHighpass = %2.1f pix   ', params.highpass ); end
            if c == 1
                preRef = preRef - imgaussfilt(preRef, params.highpass);
            end
            preData = preData - imgaussfilt(preData, params.highpass);
        end
        
        % Median filter
        if any(params.medFilter)
            if verbose,  fprintf('\nMedian filtering: [%i, %i, %i]', params.medFilter(1), params.medFilter(2), params.medFilter(3)); end
            if c == 1
                preRef = medfilt2( preRef, params.medFilter(1:2) );
            end
            preData = medfilt3( preData, params.medFilter );
        end
        %figure; subplot(1,2,1); imshow(rawRef,[]); subplot(1,2,2); imshow(preRef,[]); impixelinfo;

        % Spatial downsampling
        if params.binXY ~= 1
            if verbose, fprintf('\nSpatial downsampling factor %i', params.binXY ); end
            preData = binXY(preData, params.binXY); 
            if c == 1
                preRef = binXY(preRef, params.binXY);
            end
        end
        
        % Temporal averaging (not downsampling!)
        if params.avgT > 1 
            if verbose, fprintf('\nTemporal averaging over %i frames (sigma = %2.1f frames) ', params.avgT, params.avgTsigma); end
            [~,padChunkScans] = intersect(padScans, chunkScans); 
            avgData = cell(1,NchunkScans); % zeros(size(preData,1), size(preData,2), NchunkScans); %preData;
            parfor s = 1:NchunkScans % padChunkScans' %find(chunkScans == padScans)
                tempPadScans = padChunkScans(s)-avgTpad:padChunkScans(s)+avgTpad;
                badTempScans = find(tempPadScans< 1 | tempPadScans> NpadScans);
                tempPadScans(badTempScans) = []; 
                tempWeight = weightVector;
                tempWeight(badTempScans) = [];
                %plot(tempPadScans,  tempWeight); xlim([tempPadScans(1)-1, tempPadScans(end)+1]);
                tempWeight = repmat( permute(tempWeight, [1,3,2]), size(preData,1), size(preData,2) );
                avgData{s} = wmean(preData(:,:,tempPadScans), tempWeight, 3);  %mean(preData(:,:,tempPadScans), 3); % 
                %subplot(1,2,1); imshow( preData(:,:,padChunkScans(s)), []); 
                %subplot(1,2,2); imshow(avgData{s}, []); title(sprintf('s = %i', s));
                %pause(0.2);
            end
            avgData = cat(3, avgData{:});
            preData = avgData;
        end
        %figure; subplot(1,2,1); imshow(rawRef,[]); subplot(1,2,2); imshow(preRef,[]); impixelinfo;
        %toc
        
        % Equalize histograms to reference
        if params.histmatch
            if verbose, fprintf('\nHistogram matching enabled'); end
            %{
            figure('WindowState','maximized'); 
            subplot(2,4,1); imshow(preRef,[]); title('Reference Image (original)');
            subplot(2,4,5); histogram(preRef(:)); set(gca,'Xscale','log', 'yscale','log'); axis square
            %}
            if c == 1
                preRef = adapthisteq(uint16(preRef));
            end
            %{
            subplot(2,4,2); imshow(preRef,[]); title('Reference Image (enhanced)');
            subplot(2,4,6); histogram(preRef(:)); set(gca,'Xscale','log', 'yscale','log'); axis square
            %}
            for i = 1:NchunkScans
                tempIm = preData(:,:,i);
                %if all(matchIm == 1, 'all')
                matchIm = adapthisteq(uint16(tempIm)); %imhistmatch(uint16(tempIm), uint16(preRef));
                %else
                %    matchIm = imhistmatch(uint16(tempIm), uint16(preRef));
                %end
                preData(:,:,i) = matchIm;
                %{
                subplot(2,4,3); imshow(tempIm,[]);  title(sprintf('Original data (i = %i)', i));
                subplot(2,4,7); histogram(tempIm(:)); set(gca,'Xscale','log', 'yscale','log'); axis square;
                subplot(2,4,4); imshow(matchIm,[]);  title('Enhanced image');
                subplot(2,4,8); histogram(matchIm(:)); set(gca,'Xscale','log', 'yscale','log'); axis square;
                impixelinfo;
                pause; cla;
                %}
            end
        else
            if verbose, fprintf('\nHistogram matching disabled'); end
        end

        % Write tifs of pre-affine-registered data and reference (optional)
        if writeIntermediateTif
            % Determine tif filenames, initialize regData (optional)
            if c == 1 % only need to save reference tif once since reference is the same for all chunks
                rawRefTifPath = sprintf('%s%s_rawRef.tif', RegTempDir, saveName );
                if verbose, fprintf('\nWriting %s', rawRefTifPath); end
                WriteTiff(uint16(rawRef), rawRefTifPath); %pipe.io.writeTiff(uint16(rawRef), rawRefTifPath); %saveastiff( uint16(rawRef), rawRefTifPath, GrayOpt);
                preRefTifPath = sprintf('%s%s_preRef.tif', RegTempDir, saveName );
                if verbose, fprintf('\nWriting %s', preRefTifPath); end
                WriteTiff(uint16(preRef), preRefTifPath); %saveastiff( uint16(preRef), preRefTifPath, GrayOpt);
            end
            if Nchunk > 1
                %rawDataTifPath = sprintf('%s%s_chunk%04i_rawData.tif', RegTempDir, saveName, c  );
                preDataTifPath = sprintf('%s%s_chunk%04i_preData.tif', RegTempDir, saveName, c );
                postDataTifPath = sprintf('%s%s_chunk%04i_postData.tif', RegTempDir, saveName, c );
            else
                %rawDataTifPath = sprintf('%s%s_rawData.tif', RegTempDir, saveName );
                preDataTifPath = sprintf('%s%s_preData.tif', RegTempDir, saveName );
                postDataTifPath = sprintf('%s%s_postData.tif', RegTempDir, saveName );
            end
            % Save raw tifs
            %if verbose, fprintf('\nWriting %s', rawDataTifPath); end
            %WriteTiff(uint16(rawData), rawDataTifPath); %saveastiff( uint16(preData), preDataTifPath, GrayOpt);
            if verbose, fprintf('\nWriting %s', preDataTifPath); end
            WriteTiff(uint16(preData), preDataTifPath);
        end

        % Perform registration
        %tic
        if strcmpi(params.method, 'affine') 
            if params.turboreg
                % Use TurboReg to perform affine registration
                if isempty(alignstr)
                    [Nrow, Ncol] = size(preRef);  % Get the dimensions of the preprocessed data chunk
                    % Estimate targets the way turboreg does
                    targets = [0.5*Ncol 0.15*Nrow 0.5*Ncol 0.15*Nrow 0.15*Ncol 0.85*Nrow 0.15*Ncol 0.85*Nrow  0.85*Ncol 0.85*Nrow 0.85*Ncol 0.85*Nrow];
                    targets = round(targets);
                    targetstr = sprintf('%i ', targets);
                    % Create the text for the ImageJ macro
                    szstr = sprintf('0 0 %i %i ', Ncol - 1, Nrow - 1);
                    alignstr = sprintf('-align -window data %s -window ref %s -affine %s -hideOutput', szstr, szstr, targetstr); % -rigidBody
                end
                if ~exist('ref', 'var')
                    ref = Array2ij(preRef); %arrtoij(preRef); % 'ref' object cannot be saved
                end
                src = Array2ij(preData); %arrtoij(preData); % pipe.io.arrtoij
                if verbose, fprintf('\nRunning TurboReg... '); end
                trhl.runHL(ref, src, alignstr); % Run turboreg
                tform = trhl.getAllSourcePoints();
                targetgeotransform = trhl.getTargetPoints();
                targetgeotransform = targetgeotransform(1:3, 1:2);
                % Generate 2D transform objects for each frame
                if verbose, fprintf('\nGenerating transform objects  '); end
                for i = 1:NchunkScans % consider parallelizing here!
                    ftform = reshape(tform(i,:), 2, 3)';
                    reg_transforms{chunkScans(i)} = fitgeotrans(ftform, targetgeotransform, 'affine');
                end
            else
                % Use MATLAB's affine registration function
                tempTrans = cell(1,NchunkScans);
                parfor i = 1:NchunkScans
                    tempTrans{i} = imregtform(preData(:,:,i), preRef, params.method, optimizer, metric, 'DisplayOptimization',false, 'PyramidLevels',2); 
                    %tempTrans{i}.T
                end
                reg_transforms(chunkScans) = tempTrans;
            end
        elseif strcmpi(params.method, 'translation') || strcmpi(params.method, 'rigid') 
            % Rigid registration with MATLAB's function
            tempTrans = cell(1,NchunkScans);
            parfor i = 1:NchunkScans
                tempTrans{i} = imregtform(preData(:,:,i), preRef, params.method, optimizer, metric ); 
            end
            reg_transforms(chunkScans) = tempTrans;
        else
            fprintf('\nNo valid registration method specified');
        end
        %toc
        
        % Write out the post-registered movie to tif
        if writeIntermediateTif
            postData = zeros(size(preData)); 
            for i = 1:NchunkScans 
                postData(:,:,i) = imwarp(preData(:,:,i), reg_transforms{chunkScans(i)}, 'OutputView',imref2d(size(preRef)));
            end
            if verbose, fprintf('\nWriting %s ...', postDataTifPath); end
            WriteTiff(uint16(postData), postDataTifPath); 
        end
        
        % Correct reg_transforms for spatial downsampling
        if params.binXY ~= 1
            for i = 1:NchunkScans
                reg_transforms{chunkScans(i)}.T(3,1) = reg_transforms{chunkScans(i)}.T(3,1)*params.binXY;
                reg_transforms{chunkScans(i)}.T(3,2) = reg_transforms{chunkScans(i)}.T(3,2)*params.binXY;
            end
        end

        % Add registered frames to regData
        if nargout > 1 || writeFinalTif
            for i = 1:NchunkScans 
                unregData(:,:,K+i) = rawData(:,:,i);
                regData(:,:,K+i) = imwarp(rawData(:,:,i), reg_transforms{chunkScans(i)}, 'OutputView',imref2d(size(rawRef)));
            end
            K = K + NchunkScans;
        end
        
        % Save temporary chunk results
        chunkTforms = reg_transforms(chunkScans);
        if verbose,  fprintf('\nSaving %s', chunkSavePath); end
        if c == 1
            save(chunkSavePath, 'chunkScans','chunkTforms','params','rawRef','preRef','targets','targetstr','szstr','alignstr', '-v7.3'); % 'ref' object cannot be saved
        else
            save(chunkSavePath, 'chunkScans','chunkTforms');
        end
    else
        if verbose, fprintf('\nLoading %s', chunkSavePath); end
        load(chunkSavePath);
        reg_transforms(chunkScans) = chunkTforms;
    end
    %if Nscan > 1, waitbar(c/Nchunk, w); end
end
%if Nscan > 1, delete(w); end 
    
% Generate the final registered movie and write to tiff (optional)
if nargout > 1 || writeFinalTif
    varargout{1} = regData;
    if writeFinalTif
        unregDataTifPath = sprintf('%s%s_unregData.tif', RegTempDir, saveName );
        if verbose, fprintf('\nWriting %s ...', unregDataTifPath); end
        WriteTiff(uint16(unregData), unregDataTifPath);
        regDataTifPath = sprintf('%s%s_regData.tif', RegTempDir, saveName );
        if verbose, fprintf('\nWriting %s ...', regDataTifPath); end
        WriteTiff(uint16(regData), regDataTifPath); %toc %saveastiff( uint16(regData), regDataTifPath, GrayOpt) ;
    end
end
if verbose, fprintf('\n'); end
end