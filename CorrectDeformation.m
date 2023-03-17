raw_deform = cell(1,Nexpt);
thresh_pct = 99.9;
for x = 14 %xPresent
    % Get the results of affine registration
    [~, tformPath] = FileFinder(expt{x}.dir, 'type','mat', 'contains','regTforms');
    tformStruct = load(tformPath{1});
    regTform_raw = tformStruct.regTform;
    Nsamp = expt{x}.totScan*expt{x}.Nplane;
    tform_raw = nan(Nsamp, 8); samp = 0;
    for scan = 1:expt{x}.totScan
        for z = 1:expt{x}.Nplane
            samp = samp+1;
            tform_raw(samp,1) = scan;
            tform_raw(samp,2) = z;
            tform_raw(samp,3) = regTform_raw{z,scan}.T(3,1); % Trans AP (right/left +/-) trans_x
            tform_raw(samp,4) = regTform_raw{z,scan}.T(3,2); % Trans ML (down/up +/-)  trans_y
            tform_raw(samp,5) = regTform_raw{z,scan}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
            tform_raw(samp,6) = regTform_raw{z,scan}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
            tform_raw(samp,7) = regTform_raw{z,scan}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
            tform_raw(samp,8) = regTform_raw{z,scan}.T(2,1); % Shear ML (tilt down/right +/-) shear_y
        end
    end

    % Perform PCA, Hotelling's t-square provides a measure of how much an outlier each sample is
    [pcaCoeff, pcaScore, pcaLatent, pcaTsq, pcaExplain] = pca( tform_raw(:,3:end) ); %
    tform_raw(:,9) = pcaTsq;

    %histogram(pcaTsq, 0:0.1:100)
    Tsq = nan(expt{x}.totScan, expt{x}.Nplane);
    for samp = 1:Nsamp
        Tsq(tform_raw(samp,1), tform_raw(samp,2)) = pcaTsq(samp);
    end
    imagesc(Tsq'); caxis([0,1]); impixelinfo;

    % Determine a threshold, using the middle, pre-CSD data
    Tsq_base = Tsq(1:expt{x}.scanLims(expt{x}.preRuns(end)+1),4:expt{x}.Nplane-3); % edge planes are often sketchy
    Tsq_thresh = prctile(Tsq_base(:), thresh_pct);
    %imagesc(Tsq_base'); impixelinfo; axis image;

    figure;
    subplot(2,1,1); imagesc(Tsq');
    subplot(2,1,2); imagesc(Tsq'); caxis([0,Tsq_thresh]); impixelinfo;

    % Identify outliers
    outliers = Tsq > Tsq_thresh;
    Noutlier = sum(outliers, 2); % outliers per scan
    totOutlier = sum(Noutlier);
    outlier_frac = totOutlier/(expt{x}.totScan*expt{x}.Nplane);
    plot(Noutlier); xlim([1,expt{x}.totScan]);

    % Estimate range of inliers for each form of deformation
    inlier_mat = tform_raw(tform_raw(:,9) < Tsq_thresh, :);
    trans_range = [min(min(inlier_mat(:,3:4), [], 1)), max(max(inlier_mat(:,3:4), [], 1))];
    scale_range = [min(min(inlier_mat(:,5:6), [], 1)), max(max(inlier_mat(:,5:6), [], 1))];
    shear_range = [min(min(inlier_mat(:,7:8), [], 1)), max(max(inlier_mat(:,7:8), [], 1))];

    % Generate a corrected version of the transforms table
    tform_correct = tform_raw;
    for scan = find(Noutlier > 0)'
        % Get deformation data from that scan
        tform_scan = tform_raw(tform_raw(:,1) == scan,:);
        % Which planes need to be corrected?
        z_bad = find(outliers(scan,:));
        z_good = find(~outliers(scan,:));
        % Use interpolation to correct each form of deformation
        transAP_interp = interp1(z_good, tform_scan(z_good,3), z_bad);
        transML_interp = interp1(z_good, tform_scan(z_good,4), z_bad);
        scaleAP_interp = interp1(z_good, tform_scan(z_good,5), z_bad);
        scaleML_interp = interp1(z_good, tform_scan(z_good,6), z_bad);
        shearAP_interp = interp1(z_good, tform_scan(z_good,7), z_bad);
        shearML_interp = interp1(z_good, tform_scan(z_good,8), z_bad);

        tform_scan(z_bad,3) = transAP_interp;
        tform_scan(z_bad,4) = transML_interp;
        tform_scan(z_bad,5) = scaleAP_interp;
        tform_scan(z_bad,6) = scaleML_interp;
        tform_scan(z_bad,7) = shearAP_interp;
        tform_scan(z_bad,8) = shearML_interp;

        tform_correct(tform_raw(:,1) == scan,:) = tform_scan; % insert results into tform_correct
    end

    % How many samples were actually corrected (NaN indicates failed interpolation)
    Nfail = sum(isnan(sum(tform_correct, 2)));
    fprintf('\nFound %i outliers (%2.1f pct of samples). Corrected %i', totOutlier, outlier_frac, totOutlier-Nfail);
    nan_ind = find(isnan(tform_correct));
    tform_correct(nan_ind) = tform_raw(nan_ind);

    % Compare raw vs interpolated deformations
    transAP_raw = nan(expt{x}.totScan, expt{x}.Nplane); transAP_correct = nan(expt{x}.totScan, expt{x}.Nplane);
    transML_raw = nan(expt{x}.totScan, expt{x}.Nplane); transML_correct = nan(expt{x}.totScan, expt{x}.Nplane);
    scaleAP_raw = nan(expt{x}.totScan, expt{x}.Nplane); scaleAP_correct = nan(expt{x}.totScan, expt{x}.Nplane);
    scaleML_raw = nan(expt{x}.totScan, expt{x}.Nplane); scaleML_correct = nan(expt{x}.totScan, expt{x}.Nplane);
    shearAP_raw = nan(expt{x}.totScan, expt{x}.Nplane); shearAP_correct = nan(expt{x}.totScan, expt{x}.Nplane);
    shearML_raw = nan(expt{x}.totScan, expt{x}.Nplane); shearML_correct = nan(expt{x}.totScan, expt{x}.Nplane);
    for samp = 1:Nsamp
        scan = tform_raw(samp,1);
        z = tform_raw(samp,2);
        transAP_raw(scan,z) = tform_raw(samp,3);
        transAP_correct(scan,z) = tform_correct(samp,3);
        transML_raw(scan,z) = tform_raw(samp,4);
        transML_correct(scan,z) = tform_correct(samp,4);
        scaleAP_raw(scan,z) = tform_raw(samp,5);
        scaleAP_correct(scan,z) = tform_correct(samp,5);
        scaleML_raw(scan,z) = tform_raw(samp,6);
        scaleML_correct(scan,z) = tform_correct(samp,6);
        shearAP_raw(scan,z) = tform_raw(samp,7);
        shearAP_correct(scan,z) = tform_correct(samp,7);
        shearML_raw(scan,z) = tform_raw(samp,8);
        shearML_correct(scan,z) = tform_correct(samp,8);
    end

    figure('WindowState','maximized');
    subplot(6,2,1); imagesc(transAP_raw'); title('transAP_raw', 'Interpreter','none'); caxis(trans_range);
    subplot(6,2,2); imagesc(transAP_correct'); title('transAP_correct', 'Interpreter','none'); caxis(trans_range);
    subplot(6,2,3); imagesc(transML_raw'); title('transML_raw', 'Interpreter','none'); caxis(trans_range);
    subplot(6,2,4); imagesc(transML_correct'); title('transML_correct', 'Interpreter','none'); caxis(trans_range);
    subplot(6,2,5); imagesc(scaleAP_raw'); title('scaleAP_raw', 'Interpreter','none'); caxis(scale_range);
    subplot(6,2,6); imagesc(scaleAP_correct'); title('scaleAP_correct', 'Interpreter','none'); caxis(scale_range);
    subplot(6,2,7); imagesc(scaleML_raw'); title('scaleML_raw', 'Interpreter','none'); caxis(scale_range);
    subplot(6,2,8); imagesc(scaleML_correct'); title('scaleML_correct', 'Interpreter','none'); caxis(scale_range);
    subplot(6,2,9); imagesc(shearAP_raw'); title('shearAP_raw', 'Interpreter','none'); caxis(shear_range);
    subplot(6,2,10); imagesc(shearAP_correct'); title('shearAP_correct', 'Interpreter','none'); caxis(shear_range);
    subplot(6,2,11); imagesc(shearML_raw'); title('shearML_raw', 'Interpreter','none'); caxis(shear_range);
    subplot(6,2,12); imagesc(shearML_correct'); title('shearML_correct', 'Interpreter','none'); caxis(shear_range);
    impixelinfo;

    % Write the results back into affine2D transforms
    regTform_correct = regTform_raw;
    for samp = 1:Nsamp
        scan = tform_correct(samp,1);
        z = tform_correct(samp,2);
        regTform_correct{z,scan}.T(3,1) = tform_correct(samp,3); % Trans AP (right/left +/-) trans_x
        regTform_correct{z,scan}.T(3,2) = tform_correct(samp,4); % Trans ML (down/up +/-)  trans_y
        regTform_correct{z,scan}.T(1,1) = tform_correct(samp,5); % Scale AP (inflate/deflate >/< 1)  scale_x
        regTform_correct{z,scan}.T(2,2) = tform_correct(samp,6); % Scale ML (inflate/deflate >/< 1) scale_y
        regTform_correct{z,scan}.T(1,2) = tform_correct(samp,7); % Shear AP (tilt left/right +/-)  shear_x
        regTform_correct{z,scan}.T(2,1) = tform_correct(samp,8); % Shear ML (tilt down/right +/-) shear_y
    end

    % Save the results back into the regTforms file
    save(tformPath{1}, 'regTform_correct', 'tform_correct','tform_raw','Tsq','thresh_pct', '-append')

    % Rewrite the .sbxreg file using corrected data
    MakeSbxreg(catInfo{x});

end