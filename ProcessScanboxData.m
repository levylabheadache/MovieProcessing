%% ProcessScanboxData reads master data spreadsheet and then processes (unpacks, concatenates, registers) scanbox imaging data

% Clear any previous variables in the Workspace and Command Window to start fresh
clear; clc; close all;

% TODO --- Set the directory of where animal folders are located
dataDir = 'D:\2photon\'; %  'D:\2photon\Simone\'; %

% Parse data table

% TODO --- Set excel sheet
dataSet = 'Pollen'; % 'Macrophage'; %'AffCSD'; %  'Afferents'; % 'Macrophage'; %   'Vasculature'; %  'Astrocyte'; %  'Anatomy'; %  'Neutrophil_Simone'; %  'NGC'; % 'Neutrophil'; % 
[regParam, projParam] = DefaultProcessingParams(dataSet); % get default parameters for processing various types of data

% TODO --- Set data spreadsheet directory
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx';
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), ...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), ...
    'ref',find(contains(colNames, 'Ref')), 'edges',find(contains(colNames, 'Edge')), 'Zproj',find(contains(colNames, 'Zproj')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); loco = cell(1,Nexpt); % Tcat = cell(1,Nexpt);

% TODO --- Specify xPresent - row number(X) within excel sheet
xPresent = 7; %13; %44; %[18,22,24,30:32]; % flip(100:102); %45:47; % [66:69];
Npresent = numel(xPresent);
overwrite = false;

for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Parse data table
    [expt{x}, runInfo{x}, regParam, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir, regParam, projParam);
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x});  % , Tcat{x}

    % Get locomotion data
    for runs = expt{x}.runs % flip(expt{x}.runs) % 
        loco{x}(runs) = GetLocoData( runInfo{x}(runs), 'show',true ); 
        plot(loco{x}(regParam.refRun).Vdown)
    end

    % Use the longest period of stillness to define the reference image, or define it by hand
    % TODO -- Set the regParam.refRun and refScan to the most stable part of the pre-CSD data - IMPORTANT when you have concatenated data 
    if ~isempty(loco{x}(1).quad)
        if any(cellfun(@isempty, {loco{x}.stateDown})) %isempty(loco{x}.stateDown)
            try
                loco{x} = GetLocoState(expt{x}, loco{x}, 'dir',strcat(dataDir, expt{x}.mouse,'\'), 'name',expt{x}.mouse, 'var','velocity', 'show',true); %
            catch
                fprintf('\nGetLocoState failed for %s', expt{x}.name)
            end
        end
        % Determine reference run and scans (longest pre-CSD epoch of stillness) 
        [~,tformPath]= FileFinder(expt{x}.dir, 'contains','regTform');
        if isempty(tformPath)
            [regParam.refRun, regParam.refScan] = DetermineReference(expt{x}, Tscan{x}, loco{x}, expt{x}.preRuns, 30); % 1:4
        else
            a = load(tformPath{1}, 'params');
            regParam = a.params; clearvars a; % load previously set regParam, if it exists
        end
        % Show the scans used to define the reference
        figure;
        plot(loco{x}(regParam.refRun).Vdown); hold on; line(regParam.refScan([1,end]), [0,0], 'color','k', 'linewidth',1.5); % show the period to be used as the reference
        title(sprintf('refRun = %i', regParam.refRun)); ylabel('Velocity (cm/s)'); xlabel('Scan/Frame')
    else
        % Set reference run/scans by hand, if desired
        regParam.refRun = 1;  % plot(loco{x}(regParam.refRun).Vdown)
        regParam.refScan = 19:115; %10917:11986; % scans WITHIN the reference run
    end
    
    % Register individual runs, starting with the reference run
    if expt{x}.Nplane > 1
        runInfo{x}(regParam.refRun) = RegisterRun( runInfo{x}(regParam.refRun), regParam, 'overwrite',overwrite, 'fix',false, 'dewarp','rigid', 'interp',false);
        for runs = flip(setdiff(1:expt{x}.Nruns, regParam.refRun)) % expt{x}.runs %
            runInfo{x}(runs) = RegisterRun( runInfo{x}(runs), regParam, 'overwrite',overwrite, 'fix',false, 'dewarp','rigid', 'interp',false); %'rigid' , 'edges',[80,80,40,20]
        end
        catSbx = 'sbxdft'; % 'sbxz'; % which type of sbx to use for concatenation 
    else
        catSbx = 'sbx';  % which type of sbx files to use for concatenation
    end
    
    % Concatenate runs and metadata
    catInfo{x} = ConcatenateRunInfo(expt{x}, runInfo{x}, 'suffix','sbxcat', 'overwrite',false); % Get concatenated metadata
    ConcatenateExptRuns(expt{x}, runInfo{x}, catInfo{x}, 'sbx',catSbx); %1  expt{x}.refChan   % interRunShift = 
    catProj = WriteSbxProjection(expt{x}.sbx.cat, catInfo{x}, 'chan','green', 'type','cat', 'overwrite',overwrite, 'monochrome',true, 'RGB',true); % 
    
    % Write projections of concatenated data
    projParam.umPerPixel_target = expt{x}.umPerPixel; % avoid spatial downsampling
    projParam.edge = [60,60,20,20]; %regParam.edges; %[40,150,40,40]; %[60,60,40,20]; % [40,40,40,80];%segParams{x}.edges; %[60,60,20,20]; % [70 40 20 20]; % crop these many pixels from the [L,R,T,B] edges
    projParam.z = {22:26, 27:30}; % {3:12}; %{17:22, 24:27, 18, 19, 20, 21, 22, 23}; %{3:17, 18:28, 18, 19, 20, 21, 22, 23}; %{17:22, 23:30}; %{29:56, 5:25}; % {7:10, 18:20, 27:30};  % 1; %
    projParam.overwrite = false;
    projParam.sbx_type = {'cat'}; % , 'z'
    projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam); % write projections of unregistered data by run

    % Z interpolation (3D imaging only)
    if expt{x}.Nplane > 1
        CatInterpZ(catInfo{x}.path, catInfo{x}, 'refChan',regParam.refChan, 'refScan',regParam.refScan, 'edges',regParam.edges );
        MakeCatSbxz(catInfo{x}.path, catInfo{x}); %MakeSbxZ_new(catInfo{x}.path, catInfo{x}); % , shiftPath
        projParam.sbx_type = [projParam.sbx_type, {'z'}];
        projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam);
    end

    % Register the concatenated data (see RegisterCat3D, AlignPlanes and RegisterSBX for more info)   
    fprintf('\n   Performing planar registration... '); % (reference averaged over scans %i - %i)
    regParam.edges = projParam.edge; %  GetEdges3D( catProj, varargin )
    regParam.refScan = regParam.refScan + expt{x}.scanLims(regParam.refRun);
    %repairStruct = struct('z',1:expt{x}.Nplane, 'scan',5001:expt{x}.totScan);
    if expt{x}.Nplane > 1
        regParam = AlignPlanes( expt{x}.sbx.z, catInfo{x}, regParam, 'overwrite',overwrite, 'outPath',expt{x}.sbx.reg ); %, 'repair',repairStruct
    else
        regParam = AlignPlanes( expt{x}.sbx.cat, catInfo{x}, regParam, 'overwrite',overwrite, 'outPath',expt{x}.sbx.reg );
    end
    [~,deform{x}] = GetDeformCat3D(expt{x}, catInfo{x}, 'show',true, 'overwrite',false, 'window',find(Tscan{x}{1}<=32,1,'last'));  
    regProj = WriteSbxProjection(expt{x}.sbx.reg, catInfo{x}, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','reg', 'binT',10, 'overwrite',overwrite);
 
    % Generate downsampled, possibly z-projected, movies for each run from the concatenated data
    %[projParam.edge, x_result, y_result] = GetEdges3D( regProj(:,:,:,2), 'show',true );
    projParam.sbx_type = [projParam.sbx_type, {'reg'}];
    projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam); %  remove projParam input to load old projParam settings, if needed  

    if strcmpi(dataSet, 'afferents')
        % segment the calcium imaging data
        try
            segParams = GetSegParams( catInfo{x} );
        catch
            %zProj = projParam.z; %zProj{1} = projParam.z; % {4:10, 12:20, 22:28};
            segParams = struct('name','', 'z',[], 'edges',projParam.edge, 'seg_scan',[], 'cens_scans',1914:2099, 'min_foot',50, 'blur',2, 'corr_thresh_pct',75, 'chunk_size',30, 'hp_cutoff',51, 'min_vol',100, 'xyproj_width',100 );
            segParams.z = zProj;
        end
        segParams = SegmentSbx(catInfo{x}, segParams, 'overwrite',false); 
           
        % Generate mean/max projection images
        maxProjPath = [expt{x}.dir, expt{x}.name, '_maxProj.tif']; %[expt{x}.dir, expt{x}.name, '_affineProj.tif'];
        meanProjPath = [expt{x}.dir, expt{x}.name, '_meanProj.tif'];
        if ~exist(maxProjPath, 'file')
            [~,segProjPath] = FileFinder(expt{x}.dir, 'contains','segProj', 'type','tif'); %segProjPath = segProjPath{1};
            if ~isempty(segProjPath)
                for z = flip(1:numel(segProjPath))
                    tempSegMov{z} = imread_big(segProjPath{z});
                end
                tempSegMov = cat(4, tempSegMov{:});
                expt{x}.maxProj = max(tempSegMov, [] , [3,4]);
                saveastiff(expt{x}.maxProj, maxProjPath);
                expt{x}.meanProj = mean(tempSegMov, [3,4]);
                saveastiff(uint16(expt{x}.meanProj), meanProjPath);
                clearvars tempSegMov;
            else
                fprintf('\nseg projections do not exist - cannot generate mean/max projections!');
            end
        else
            expt{x}.maxProj = loadtiff( maxProjPath );
            expt{x}.meanProj = loadtiff( meanProjPath );
        end
        
        % Make the final ROIs from the segmentation data
        ROI{x} = MakeROI3D(expt{x}, 'overwrite',false, 'corrPrct',75, 'minFoot',50);  % , 'corrPrct',90, [ROI{x}, preROI{x}] overwriteROI
        expt{x}.Nroi = numel(ROI{x}); % Nauto(x)
        [expt{x}.roiProj, ~, expt{x}.roiLabel] = VisualizeSegmentation(expt{x}, ROI{x}, 'overwrite',false);   %  
        % Extract/process fluor signals
        for runs = flip(expt{x}.runs)
            fluor{x}(runs) = GetFluor3D(runInfo{x}(runs), 'overwrite',false);  % expt{x}, r
        end
        fluor{x} = GetROIfluor(expt{x}, catInfo{x}, ROI{x}, fluor{x}, deform{x}, loco{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'lp',0, 'deconvolve',false, 'overwrite',false); 
        % Show the final results
        defVars = {'transAP', 'transML', 'transMag', 'scaleAP', 'scaleML', 'scaleMag', 'stretchAP', 'stretchML', 'shearAP', 'shearML', 'shearMag', 'shiftZ', 'DshiftZ'}; %
        NdefVars = numel( defVars ); %, 'dShiftZ'
        allVars = [{'fluor'}, defVars, {'velocity'}]; NallVars = numel(allVars);
        viewLims = struct('trans',[-Inf,Inf], 'scale',[-Inf,Inf], 'stretch',[-Inf,Inf], 'shear',[-0.03,0.03], 'shift',[-3, 3], 'velocity',[-3, 15], 'fluor',[-5, 5]); 
        ViewResults3D( expt{x}, Tscan{x}, deform{x}, loco{x}, fluor{x}, allVars, ROI{x}, 'cat',true, 'limits',viewLims); %
    end
end
clearvars Expt;