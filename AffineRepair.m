sbxInfo = catInfo(x);
[deformCat, deform{x}, affParams, badInd] = GetDeformCat3D( catInfo(x), deformLim, 'show',true, 'overwrite',true);  %   [deformCat, deform{x}, affParams, badInd] = 

%plot(deformCat.RS_final, deformCat.CS_final, '.');

% Load ref and current results
sbxInputPath = sbxInfo.path; %sbxPath{runs};
[fDir, fName] = fileparts(sbxInputPath);
refTifPath = sprintf('%s%s%s_affineRef.tif', fDir, '\', fName);
refVol = loadtiff(refTifPath);
tformPath = strcat(fDir,'\',fName,'_affine_tforms.mat');
fprintf('\nLoading %s... \n', tformPath );
load( tformPath, 'tforms_all', 'params'); %  , 'refVol'   'sbxPath', 'sbxInfo', 'affParams', 'tforms_all'

%params = affParams;
tifDir = strcat(fDir,'\AffTemp\'); mkdir(tifDir);

%% 

csdScans = csdBout{x}(expt(x).csd).scan{1} + expt(x).scanLims(expt(x).csd);
WriteSbxPlaneTif(expt(x).sbx, catInfo(x), 1, 'firstScan',csdScans(1), 'Nscan',numel(csdScans), 'dir','D:\2photon\DL72\170614_DL72\AffTemp\', 'name','DL72_CSD', 'verbose',true, 'overwrite',true);
WriteSbxPlaneTif(expt(x).sbx, catInfo(x), 1, 'firstScan',57280, 'Nscan',1400, 'dir','D:\2photon\DL72\170614_DL72\AffTemp\', 'name','DL72_quiescent', 'verbose',true, 'overwrite',true);
%% Change affine parameters and  repair
params.lowpass = 0;
params.highpass = 0;
params.medFilter = [9,9,1]; % ,3
params.histmatch = false; % true; %%true;
%params.medFilter = [0,0,0]; %[3,3,3]; 

% Repair bad registration
tforms_repair = tforms_all;
zRepair = 14; % % 1; %[16,27:29];
for z = zRepair %zRepair
    %tempStretchML = vertcat(deform{x}.stretchML);
    repairScan = find(abs(deformCat.scaleAP(:,z)) > 1.01 | abs(deformCat.scaleAP(:,z)) < 0.91)'; % 1:1860; %1860:1959; %find(abs(deformCat.stretchAP(:,zRepair)) > 50 )'; %find(deformCat.shearAP > 0.01)'; % find(abs(tempStretchML(:,zRepair)) > 80 )'; %  %350; 1893:1894; %
    %repairScans = csdBout{x}(run).scan{1} + expt(x).scanLims(expt(x).csd);
    NrepairScan = numel(repairScan);
    c = 0;
    w = waitbar(0, sprintf('Repairing plane %i', z) );
    for s = repairScan 
        c = c + 1;
        try
            temp_tform = AffineTurboReg(sbxInputPath, sbxInfo, refVol(:,:,z), params, z, 'firstScan',s, 'Nscan',1, ...
                'dir',tifDir, 'name',sprintf('%s_z%02d', fName, z), 'verbose',false, 'intTif',false, 'finalTif',false, 'overwrite',true); 
            
            
            %temp_tform = AffineTurboReg(sbxInputPath, sbxInfo, refVol(:,:,z), params, z, 'firstScan',1, 'Nscan',300, ...
            %    'dir',tifDir, 'name',sprintf('%s_z%02d', fName, z), 'verbose',false, 'intTif',false, 'finalTif',true, 'overwrite',true); 
            
            
            tforms_repair(z,s) = temp_tform(1,s); % (z,s)
            waitbar(c/NrepairScan, w); % 
        catch
            fprintf('\nRepair failed: [z,s] = [%i, %i]', z, s);
        end
    end
    delete(w);
    
    for s = 1:sbxInfo.totScan
        deformOrig.transAP(s,z) = tforms_all{z,s}.T(3,1); % Trans AP (right/left +/-) trans_x
        deformOrig.transML(s,z) = tforms_all{z,s}.T(3,2); % Trans ML (down/up +/-)  trans_y
        deformOrig.scaleAP(s,z) = tforms_all{z,s}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
        deformOrig.scaleML(s,z) = tforms_all{z,s}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
        deformOrig.shearAP(s,z) = tforms_all{z,s}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
        deformOrig.shearML(s,z) = tforms_all{z,s}.T(2,1); % Shear ML (tilt down/right +/-) shear_y

        deformRepair.transAP(s,z) = tforms_repair{z,s}.T(3,1); % Trans AP (right/left +/-) trans_x
        deformRepair.transML(s,z) = tforms_repair{z,s}.T(3,2); % Trans ML (down/up +/-)  trans_y
        deformRepair.scaleAP(s,z) = tforms_repair{z,s}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
        deformRepair.scaleML(s,z) = tforms_repair{z,s}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
        deformRepair.shearAP(s,z) = tforms_repair{z,s}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
        deformRepair.shearML(s,z) = tforms_repair{z,s}.T(2,1); % Shear ML (tilt down/right +/-) shear_y
    end
end
%%
for z = zRepair
    % {

    %}
    close all;
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    sp(1) = subplot(2,2,1);
    plot( deformOrig.scaleAP(:,z) ); % hold on;
    title( sprintf('z = %i:  Original',z));
    sp(2) = subplot(2,2,3);
    plot( deformRepair.scaleAP(:,z));
    title('Repair');
    ylim([0.95,1.05]); xlim([repairScan(1)-5, repairScan(end)+5]);
    subplot(2,2,[2,4]); % sp(3) = 
    plot( deformOrig.scaleML(repairScan,z), deformRepair.scaleML(repairScan,z), '.' );
    xlabel('Original'); ylabel('Repair')
    axis square;
    linkaxes(sp,'xy')
    pause;
end



%%
tforms_all = tforms_repair;
save(tformPath, 'sbxInputPath', 'sbxInfo', 'refVol', 'params', 'tforms_all', '-mat');

%}


%% Repair by interpolation

for z = zRepair
    tempTransAP = deformOrig.transAP(:,z); tempTransAP(repairScan) = NaN;
    tempTransML = deformOrig.transML(:,z); tempTransML(repairScan) = NaN;
    tempScaleAP = deformOrig.scaleAP(:,z); tempScaleAP(repairScan) = NaN;
    tempScaleML = deformOrig.scaleML(:,z); tempScaleML(repairScan) = NaN;
    tempShearAP = deformOrig.shearAP(:,z); tempShearAP(repairScan) = NaN;
    tempShearML = deformOrig.shearML(:,z); tempShearML(repairScan) = NaN;
    
    goodScan = 1:size(tempTransAP,1);
    goodScan(repairScan) = [];
    
    tempTransAP(repairScan) = interp1( goodScan, tempTransAP(goodScan), repairScan, 'spline' );
    tempTransML(repairScan) = interp1( goodScan, tempTransML(goodScan), repairScan, 'spline' );
    tempScaleAP(repairScan) = interp1( goodScan, tempScaleAP(goodScan), repairScan, 'spline' );
    tempScaleML(repairScan) = interp1( goodScan, tempScaleML(goodScan), repairScan, 'spline' );
    tempShearAP(repairScan) = interp1( goodScan, tempShearAP(goodScan), repairScan, 'spline' );
    tempShearML(repairScan) = interp1( goodScan, tempShearML(goodScan), repairScan, 'spline' );
    
    for s = repairScan
        tforms_repair{z,s}.T(3,1) = tempTransAP(s);
        tforms_repair{z,s}.T(3,2) = tempTransML(s);
        tforms_repair{z,s}.T(1,1) = tempScaleAP(s);
        tforms_repair{z,s}.T(2,2) = tempScaleML(s);
        tforms_repair{z,s}.T(1,2) = tempShearAP(s);
        tforms_repair{z,s}.T(2,1) = tempShearML(s);
    end
end
tforms_all = tforms_repair;
save(tformPath, 'sbxInputPath', 'sbxInfo', 'refVol', 'params', 'tforms_all', '-mat');

%% Apply the transforms and generate sbx_affine

Nplane = sbxInfo.Nplane; 
Nscan = sbxInfo.totScan; 
Nrow = sbxInfo.sz(1); 
Ncol = sbxInfo.sz(2); 
Nchan = sbxInfo.nchan;

sbxAffPath = strcat( fDir, '\', fName, '.sbx_affine' );

tic
w = waitbar(0,'writing .sbx\_affine');
rw = pipe.io.RegWriter(sbxAffPath, sbxInfo, '.sbx_affine', true);
if Nchan == 1
    fprintf('\n     Writing monochrome affine-registered sbx file'); 
    aff_chunk = zeros(Nrow, Ncol, Nplane);
    for s = 1:Nscan
        % Apply affine transformation to each plane of the scan
        input_chunk = readSBX(sbxInputPath, sbxInfo, Nplane*(s-1)+1, Nplane, params.refChan, []); % readSBX(sbxInputPath, sbxInfo, s, 1, params.refChan, []); %pipe.imread(sbxInputPath, Nplane*(s-1)+1, Nplane, params.refChan, []);
        for z = 1:Nplane % parfor is actually slower
            aff_chunk(:,:,z) = imwarp(input_chunk(:,:,z), tforms_repair{z,s}, 'OutputView',imref2d([Nrow,Ncol]));
        end
        % Write the results to sbx_affine
        rw.write(squeeze(uint16(aff_chunk))); %rw.write(squeeze(uint16(tempScan)));
        waitbar(s/Nscan, w);
    end
else
    fprintf('\n     Writing two-color affine-registered sbx file'); 
    for s = 1:Nscan
        input_chunk = pipe.imread(sbxInputPath, Nplane*(s-1)+1, Nplane, -1, []); %raw_scan = pipe.imread(path, Nplane*(s-1)+1, Nplane, -1, []);
        input_chunk = input_chunk(:,params.edges(3)+1:end-params.edges(4), params.edges(1)+1:end-params.edges(2),:);
        input_chunk = reshape(input_chunk, Nchan, NxCrop, NyCrop, Nplane, []);
        trans_scan = zeros( Nchan, NxCrop, NyCrop, Nplane );
        for c = 1:2
            for z = flip(1:Nplane)
                trans_scan(c,:,:,z) = imwarp( squeeze(input_chunk(c,:,:,z)), tforms_repair{z,s}, 'OutputView',imref2d([Nrow,Ncol])); 
            end
        end
        % Write to sbx_affine
        tempScan = zeros(2, Nrow, Ncol, Nplane, 'uint16');
        tempScan(:,params.edges(3)+1:end-params.edges(4), params.edges(1)+1:end-params.edges(2),:) = trans_scan;
        rw.write( tempScan ); % rw.write(squeeze(uint16(tempScan)));
        waitbar(s/Nscan, w);
    end
end
rw.delete;
delete(w);
toc
%
ZtifDir = strcat(sbxInfo.dir, 'Ztifs', '\'); mkdir( ZtifDir );
affineStack = WriteSbxPlaneTif(sbxAffPath, sbxInfo, 1, 'dir',ZtifDir, 'name',fName, 'type','aff', 'binT',10, 'verbose',true, 'chan',1, 'overwrite',true ); % , 'edge',params.edges, 'scale',1
affineMean = mean(affineStack, 3);
affProjPath = strcat(fDir, '\', fName, '_affineProj.tif' ); % strcat(pathTemplate,'_affineProj.tif');
saveastiff(uint16(affineMean), affProjPath );
