function refShift = CorrectData2D(sbxPath, sbxInfo, regParams, shiftPath, sbxDftPath, varargin)
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'sbxInfo', @isstruct )
addRequired(IP, 'regParams', @isstruct )
addRequired(IP, 'shiftPath', @ischar )
addRequired(IP, 'sbxDftPath', @ischar )
addParameter(IP, 'sigma', [0,0,0], @isnumeric)
addParameter(IP, 'verbose', true, @islogical)
addParameter(IP, 'overwrite', false, @islogical)
parse(IP, sbxPath, sbxInfo, regParams, shiftPath, sbxDftPath, varargin{:}); % tforms_optotune,
sigma3 = IP.Results.sigma;
verbose = IP.Results.verbose;
overwrite = IP.Results.overwrite;
if ~exist(sbxDftPath,'file') || overwrite
    if ~exist(shiftPath, 'file') || overwrite
        % Load the cropped data and calculate Fourier transforms of each scan
        cropData = WriteSbxPlaneTif(sbxPath, sbxInfo, 1, 'chan',regParams.refChan, 'edges',regParams.edges, 'monochrome',false, 'RGB',false, 'overwrite',true, 'verbose',verbose);  %overwrite
        % , 'firstScan',460 , 'Nscan',1400 , 'scale',regParams.binXY
        %WriteTiff(cropData, sprintf('%s%s_crop.tif', sbxInfo.dir, sbxInfo.exptName));
        if any(sigma3 ~= 0) 
            if verbose, fprintf('\nApplying gaussian filtering, sigma [x,y,z] = [%i, %i, %i] pix or frames', sigma3);  end
            cropData = imgaussfilt3(cropData, sigma3, 'padding','symmetric' );
            %WriteTiff(cropData, sprintf('%s%s_sigma%i_%i_%i.tif', sbxInfo.dir, sbxInfo.exptName, sigma3));
        end

        % Calculate FFTs of each frame
        if verbose, fprintf('\nCalculating FFTs'); end
        dataSize = size(cropData);
        cropFT = zeros(dataSize);
        for s = 1:sbxInfo.Nscan % par is slower
            cropFT(:,:,s) = fft2(cropData(:,:,s));
        end

        % Find the window of scans with minimal motion and use that window to generate a reference image
        if verbose, fprintf('\nGenerating the reference image'); end
        scanShift = zeros(size(cropData,3)-1, 4); % shift relative to prior scan
        for s = 1:dataSize(3)-1
            scanShift(s,:) = dftregistration(cropFT(:,:,s),  cropFT(:,:,s+1), 10);
        end
        scanShiftMag = sum(scanShift(:,3:4).^2, 2); 
        %figure; plot( scanShiftMag );
        scanShiftMean = movmean(scanShiftMag, regParams.window); 
        minShiftScan = find(scanShiftMean == min(scanShiftMean));
        refScan = median(minShiftScan); 
        %figure;  plot( scanShiftMean ); hold on; plot(minShiftScan, scanShiftMean(minShiftScan), 'ro'); plot(refScan, scanShiftMean(refScan), 'ko')
        % Construct the reference image from a window around refScan
        refWindow = round(refScan-regParams.window/2):round(refScan+regParams.window/2);
        refWindow(refWindow < 1 | refWindow > sbxInfo.Nscan) = [];
        refImage = mean(cropData(:,:,refWindow), 3); %imshow(refImage, [])
        clearvars cropData;
        refFT = fft2(refImage);

        % Register data to reference image
        if verbose, fprintf('\nPerforming DFT registration'); end
        %tic
        refShift = zeros(dataSize(3), 4); % shift relative to prior scan
        parfor s = 1:dataSize(3) % par is faster
            refShift(s,:) = dftregistration(refFT,  cropFT(:,:,s), 10);
        end
        %toc
        refShift = refShift(:,[4,3]); % for compatibility with imtranslate  regParams.binXY*
        refShiftMag = sqrt(sum(refShift.^2, 2)); % plot( refShiftMag );
        refShiftMean = movmean(refShiftMag, regParams.window);
        %figure; plot( refShiftMean );
        fprintf('\nSaving %s', shiftPath)
        save(shiftPath, 'refShift','refShiftMag', 'refShiftMean', '-mat', '-v7.3');
    else
        fprintf('\nLoading %s', shiftPath)
        load(shiftPath);
    end

    % Write the sbxdft file
    [chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.Nscan, 'N',10);
    fprintf('\n     Writing %s\n', sbxDftPath); tic
    rw = SbxWriter(sbxDftPath, sbxInfo, '.sbxdft', true); % pipe.io.RegWriter(catSbxPath, catInfo, catExt, true);
    w = waitbar(0, 'writing .sbxdft');
    if verbose, tic; end
    if sbxInfo.nchan == 2 % Single plane, multi color
        for chunk = 1:Nchunk
            data_chunk = readSBX(sbxPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), -1, []);
            data_chunk = permute(data_chunk, [2,3,4,1]);
            for s = 1:chunkLength(chunk)
                for chan = 1:2
                    data_chunk(:,:,s,chan) = imtranslate( data_chunk(:,:,s,chan), refShift(chunkLims(chunk,1)+s-1,:)  ); %imwarp(data_chunk(:,:,s,chan), regTform{1,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                end
            end
            data_chunk = permute(data_chunk, [4,1,2,3]);
            rw.write(data_chunk);
            waitbar(chunk/Nchunk, w);
        end
    else % Single plane, single color
        usePMTind = DeterminePMT('both', sbxInfo);
        for chunk = 1:Nchunk
            data_chunk = readSBX(sbxPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), usePMTind, []);
            for s = 1:chunkLength(chunk)
                data_chunk(:,:,s) = imtranslate( data_chunk(:,:,s), refShift(chunkLims(chunk,1)+s-1,:) );
            end
            rw.write(data_chunk);
            waitbar(chunk/Nchunk, w);
        end
    end
    rw.delete;
    delete(w);
    if verbose, tic; end
else
    fprintf('\n%s already exists and overwrite is disabled - returning\n', sbxDftPath);
    % GITHUB TEST
end
end