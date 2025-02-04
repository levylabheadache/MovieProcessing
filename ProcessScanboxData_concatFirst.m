%% ProcessScanboxData reads master data spreadsheet and then processes (unpacks, concatenates, registers) scanbox imaging data

% Clear any previous variables in the Workspace and Command Window to start fresh
clear; clc; close all;

% TODO --- Set the directory of where animal folders are located
dataDir =  'D:\2photon\Simone\Simone_Vasculature\'; %  'D:\2photon\Simone\Simone_Macrophages\'; %  'Z:\2photon\Simone\Simone_Macrophages\'; %  'D:\2photon\'; %  'D:\2photon\Simone\'; % 

% Parse data table

% TODO --- Set excel sheet
dataSet = 'Vasculature'; %'Macrophage'; % 'AffCSD'; %  'Pollen'; 'Vasculature'; 'Astrocyte'; %  'Anatomy'; %  'Neutrophil_Simone'; %  'NGC'; % 'Neutrophil'; % 'Afferents'
[regParam, projParam] = DefaultProcessingParams(dataSet); % get default parameters for processing various types of data

regParam.method = 'rigid';
regParam.name = 'rigid';

% TODO --- Set data spreadsheet directory
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx';
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), 'vascChan',find(contains(colNames, 'VascChan')),...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), ...
    'ref',find(contains(colNames, 'Ref')), 'edges',find(contains(colNames, 'Edge')), 'Zproj',find(contains(colNames, 'Zproj')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); loco = cell(1,Nexpt); % Tcat = cell(1,Nexpt);

% TODO --- Specify xPresent - row number(X) within excel sheet
xPresent = 70; %; %13; %44; %[18,22,24,30:32]; % flip(100:102); %45:47; % [66:69];
Npresent = numel(xPresent);
overwrite = false;

for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Parse data table
    [expt{x}, runInfo{x}, regParam, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir, regParam, projParam);
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x});  % , Tcat{x}

    % Get locomotion data
    for runs = expt{x}.runs % flip(expt{x}.runs) % 
        loco{x}(runs) = GetLocoData( runInfo{x}(runs), 'show',true ); 
        %plot(loco{x}(regParam.refRun).Vdown)
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
        refScanFig = figure;
        plot(loco{x}(regParam.refRun).Vdown); hold on; line(regParam.refScan([1,end]), [0,0], 'color','k', 'linewidth',1.5); % show the period to be used as the reference
        title(sprintf('refRun = %i', regParam.refRun)); ylabel('Velocity (cm/s)'); xlabel('Scan/Frame')
        savefig(refScanFig, fullfile(expt{x}.dir,['refScanFig.fig']));

        regParam.refScan = regParam.refScan + expt{x}.scanLims(regParam.refRun);
    else
        % Set reference run/scans by hand, if desired
        %regParam.refRun = 1;  % plot(loco{x}(regParam.refRun).Vdown)
        %regParam.refScan = 19:115; %10917:11986; % scans WITHIN the reference run
        regParam.refScan = 90:400;
    end
    
    % Concatenate unprocessed runs and metadata
    catInfo{x} = ConcatenateRunInfo(expt{x}, runInfo{x}, 'suffix','sbxcat', 'overwrite',false); % Get concatenated metadata
    if expt{x}.Nplane > 1
        ConcatenateRuns(expt{x}, runInfo{x}, catInfo{x}, 'sbx','sbxfix'); %1  expt{x}.refChan   % interRunShift = 
    else
        ConcatenateRuns(expt{x}, runInfo{x}, catInfo{x}, 'sbx','sbx'); 
    end
    catProj = WriteSbxProjection(expt{x}.sbx.cat, catInfo{x}, 'chan','both', 'type','cat', 'overwrite',overwrite, 'monochrome',true, 'RGB',true); % 
    projParam.sbx_type = {'cat'}; % , 'z'

    % Correct optotune warping
    if expt{x}.Nplane > 1 && ~exist(expt{x}.sbx.opt, 'file')
        DewarpSbx(expt{x}, catInfo{x}, regParam); % dewarpType
        WriteSbxProjection(expt{x}.sbx.opt, catInfo{x}, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','opt', 'overwrite',overwrite);
        projParam.sbx_type = [projParam.sbx_type, {'opt'}]; % 
    end

    % Calculate rigid corrections
    CorrectExpt(expt{x}, catInfo{x}, regParam);
    WriteSbxProjection(expt{x}.sbx.dft, catInfo{x}, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','dft', 'overwrite',overwrite);
    projParam.sbx_type = [projParam.sbx_type, {'dft'}]; % 

    % Z interpolation (3D imaging only)
    if expt{x}.Nplane > 1
        ExptInterpZ(expt{x}, catInfo{x}, regParam ); % .sbx.dft
        %MakeCatSbxz(catInfo{x}.path, catInfo{x}); %MakeSbxZ_new(catInfo{x}.path, catInfo{x}); % , shiftPath
        WriteSbxProjection(expt{x}.sbx.z, catInfo{x}, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','z', 'overwrite',overwrite);
        projParam.sbx_type = [projParam.sbx_type, {'z'}];
    end

    % Write projections of concatenated data at various stages of processing
    %TODO: set the projParam.z to the frames that you would like to process
    projParam.umPerPixel_target = expt{x}.umPerPixel; % set to expt{x}.umPerPixel to avoid spatial downsampling, otherwise give a number
    projParam.edge = [80,80,40,40];  % [80,80,40,40], crop these many pixels from the [L,R,T,B] edges
    %projParam.z = {5, 5:7}; %{4:5, 5:7}{5, 7:8} {3:12}; %{17:22, 24:27, 18, 19, 20, 21, 22, 23}; %{3:17, 18:28, 18, 19, 20, 21, 22, 23}; %{17:22, 23:30}; %{29:56, 5:25}; % {7:10, 18:20, 27:30};  % 1; %
    projParam.overwrite = false;
    %projParam.sbx_type = {'cat','dft','z','reg'}; % , 'z'
    projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam); % write projections of unregistered data by run

    % Register the concatenated data (see RegisterCat3D, AlignPlanes and RegisterSBX for more info)   
    fprintf('\n   Performing planar registration... '); % (reference averaged over scans %i - %i)
    %regParam.refScan = regParam.refScan + expt{x}.scanLims(regParam.refRun);
    %repairStruct = struct('z',1:expt{x}.Nplane, 'scan',5001:expt{x}.totScan);
    %
    if expt{x}.Nplane > 1
        regParam = AlignPlanes( expt{x}.sbx.z, catInfo{x}, regParam, 'overwrite',overwrite, 'outPath',expt{x}.sbx.reg ); %, 'repair',repairStruct
    else
        regParam = AlignPlanes( expt{x}.sbx.dft, catInfo{x}, regParam, 'overwrite',overwrite, 'outPath',expt{x}.sbx.reg );
    end

    %Deformation Limits
    deformLim = struct('trans',[-Inf,Inf], 'scale',[0.95, 1.05], 'shear',[-0.03, 0.03], 'shift',[-3.5, 3.5]); %, 'stretch' ,100*[-1,1]); %'trans',[-0.5, 0.5], scale(%), trans(px), z(planes), shear(
    [~,deform{x}, ~, badInd] = GetDeformCat3D(expt{x}, catInfo{x}, deformLim, 'show',true, 'overwrite',true, 'window',find(Tscan{x}{1}<=32,1,'last'));  
    regProj = WriteSbxProjection(expt{x}.sbx.reg, catInfo{x}, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'type','reg', 'overwrite',overwrite); % , 'binT',10, 'overwrite',overwrite
    projParam.sbx_type = [projParam.sbx_type, {'reg'}];

    % Generate downsampled, possibly z-projected, movies for each run from the concatenated data
    projParam = GenerateExptProjections(expt{x}, catInfo{x}, Tscan{x}, projParam); %  remove projParam input to load old projParam settings, if needed     

    % Locomotion bouts
    for runs = 1:expt{x}.Nruns
        [periBout{x}(runs), periParam{x}(runs), ~] = PeriLoco(expt{x}, Tscan{x}{runs}, loco{x}(runs), deform{x}(runs), 'base',10, 'min_vel_on',2, 'merge',true,'show',true); % 'iso',[10,0] , showDefVars , fluor{x}(runs).F.axon, 'run',2,  
    end
end
clearvars Expt;