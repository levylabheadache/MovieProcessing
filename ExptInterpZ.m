function ExptInterpZ(expt, sbxInfo, regParam, varargin)
IP = inputParser;
addRequired(IP, 'expt', @isstruct )
addRequired(IP, 'sbxInfo', @isstruct )
addRequired(IP, 'regParam', @isstruct )
addParameter(IP, 'chunkSize', 15, @isnumeric )
addParameter(IP, 'refType', 'mean', @ischar); %options are 'median' or 'mean' for projecting the reference volume
addParameter(IP, 'planescorr', 3, @isnumeric); %how many planes up and down to search for Z DFT reg
addParameter(IP, 'blur', 1, @isnumeric); %width of gaussian to blur for DFT reg
addParameter(IP, 'keep',0.95, @isnumeric); %what proportion of frame to accountf or with shifts
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, sbxInfo, regParam, varargin{:} );
[usePMT, ~] = DeterminePMT(regParam.refChan, sbxInfo);
planes_corr = IP.Results.planescorr;
overwrite = IP.Results.overwrite;
% save the results
zinterpPath = strcat(sbxInfo.dir, sbxInfo.exptName,'_zinterp.mat'); %IP.Results.shiftpath;
if ~exist(zinterpPath,'file') || overwrite
    %Ncheck = numel(planescorr+1:sbxInfo.Nplane-planescorr);
    % Generate reference volume
    fprintf('\nAveraging scans %i - %i for reference volume\n', regParam.refScan(1), regParam.refScan(end));
    ref_vol = WriteSbxProjection(expt.sbx.dft, sbxInfo, 'firstScan',regParam.refScan(1), 'Nscan',numel(regParam.refScan), 'type','zref', 'chan',regParam.refChan, 'verbose',true, 'monochrome',true, 'overwrite',false); % , 'regParam.edges',regParam.edges , 'regParam.binXY',regParam.binXY
    ref_vol = ref_vol(regParam.edges(3)+1:end-regParam.edges(4),regParam.edges(1)+1:end-regParam.edges(2),:,:);
    ref_vol = imresize(ref_vol,1/regParam.binXY);
    corr_length = 2*planes_corr+1;
    % Interpolate each scan to the reference
    fprintf('\nEstimating z-shift relative to reference volume')
    w = waitbar(0,'Estimating z-shift relative to reference volume...');
    RS = zeros(sbxInfo.Nplane, sbxInfo.totScan); CS = zeros(sbxInfo.Nplane, sbxInfo.totScan); ZS = zeros(sbxInfo.Nplane, sbxInfo.totScan); % row shift, column shifts, z shifts
    for scan = 1:sbxInfo.totScan
        % Load the current scan
        %curr_scan = WriteSbxProjection(expt.sbx.dft, sbxInfo, 'firstScan',scan, 'Nscan',1, 'chan',regParam.refChan, 'verbose',false, 'monochrome',false, 'regParam.edges',regParam.edges, 'regParam.binXY',regParam.binXY, 'overwrite',false); %
        curr_scan = readSBX(expt.sbx.dft, sbxInfo, scan, 1, usePMT); % (path, info, firstScan, Nscan, pmt, z)
        curr_scan = curr_scan(regParam.edges(3)+1:end-regParam.edges(4),regParam.edges(1)+1:end-regParam.edges(2),:,:);
        curr_scan = imresize(curr_scan,1/regParam.binXY);
        for z = planes_corr+1:sbxInfo.Nplane-planes_corr % j = considered plane
            corr_z = nan(1,corr_length);
            for i = z-planes_corr:z+planes_corr % i = reference plane
                corr_z(1,i-z+planes_corr+1) = corr2(ref_vol(:,:,i), curr_scan(:,:,z));
            end
            [~,J] = max(corr_z);
            % Set interpolation vectors and degree regarding matrix size
            if  J-5 > 0 && J+5 <= corr_length
                idx = 5;
            else
                idx = min([corr_length-J, J-1]);
            end
            % Find Z shift based on spatial correlations
            x = J-idx:0.01:J+idx;
            FitOrder = idx;
            P = polyfit(J-idx:J+idx, corr_z(J-idx:J+idx),FitOrder);
            CorrelationFit = polyval(P, x);
            [~,x_max_ind] = max(CorrelationFit); % max of the polynomial curve
            % recalcutate x and y shifts
            output = dftregistrationAlex(fft2(ref_vol(:,:,z)), fft2(curr_scan(:,:,z)), 100);
            row_shift(z,1) = output(1);
            column_shift(z,1) = output(2);
            z_shift(z,1) = x(x_max_ind)-planes_corr-1;
        end
        z_shift_length = numel(z_shift);
        % add shifts into output matrix
        RS(:,scan) = cat(1,row_shift, zeros(planes_corr,1));
        CS(:,scan) = cat(1,column_shift, zeros(planes_corr,1));

        % ensuring strict monotony, necessary?
        psteps = (1:z_shift_length)'; %ones(z_shift_length,1).*(1:z_shift_length)';
        zaux = -z_shift + psteps;
        count = 0;
        while ~issorted(zaux) && count < 10
            count = count + 1;
            for plane = 2:z_shift_length
                if zaux(plane) - zaux(plane-1) <= 0
                    if ismember(count, [1,2,5,6,9,10]) % mod(count,2)==1;
                        zaux(plane-1) = NaN;
                    else
                        zaux(plane) = NaN;
                    end
                end
            end
            zaux = naninterp(zaux);
        end
        z_shift = zaux - psteps;
        ZS(:,scan) = cat(1,-z_shift, zeros(planes_corr,1));
        waitbar(scan/sbxInfo.totScan, w)
    end
    delete(w);
    fprintf('\nSaving %s', zinterpPath)
    save(zinterpPath, 'RS','CS','ZS','regParam','expt','sbxInfo', '-mat'); % 'RS1','CS1','ZS1','RS2','CS2','RS3','CS3','RS_chunk','CS_chunk','ZS_chunk','regParam.binXY'
    
else
    fprintf('\nLoading %s', zinterpPath); % fprintf('\n%s already exists\n', shiftPath);
    load(zinterpPath)
end

% Check if the sbxz file already exists
[~, sbxzPath] = FileFinder( sbxInfo.dir, 'type','mat', 'type','sbxz' );
if isempty(sbxzPath)
    % Load interpolation results
    [~, zinterpPath] = FileFinder( sbxInfo.dir, 'type','mat', 'contains','zinterp' );
    if ~isempty(zinterpPath)
        fprintf('\nLoading %s', zinterpPath{1})
        load(zinterpPath{1},'-mat');
        % Generate .sbxz file
        rw = SbxWriter(expt.sbx.z, sbxInfo, '.sbxz');
        w = waitbar(0,'Generating sbxz');
        for scan = 1:sbxInfo.totScan
            % Load one scan at a time
            data_scan = readSBX(expt.sbx.dft, sbxInfo, scan, 1, -1, []);
            % Apply Zshifts
            if sbxInfo.nchan == 1
                data_scan = ApplyZShiftInterpolateFBS(data_scan, ZS(:,scan), CS(:,scan), RS(:,scan)); %toc % apply shifts, this step is slow
                data_scan = reshape(data_scan, [size(data_scan,[1,2]), prod(size(data_scan,[3,4]))] );
            else
                data_scan = permute(data_scan, [2,3,4,5,1]); % ApplyZShiftInterpolateFBS expects spatial dimensions first
                for c = 1:2 % (parallel seems to be slower than serial)
                    data_scan(:,:,:,:,c) = ApplyZShiftInterpolateFBS(data_scan(:,:,:,:,c), ZS(:,scan), CS(:,scan), RS(:,scan));
                end
                data_scan = permute(data_scan, [5,1,2,3,4]); %toc % regwriter expects color dim first
                data_scan = reshape(data_scan, [size(data_scan,[1,2,3]), prod(size(data_scan,[4,5]))] );
            end
            % Write result to .sbxz
            rw.write(data_scan);
            waitbar(scan/sbxInfo.totScan, w);
        end
        toc
        delete(w);
        rw.delete;
    else
        error('\nNo zinterp file found!');
    end
else
    fprintf('\n%s already exists - returning\n', sbxzPath{1})
end