function [deformResultStruct, deformSubStruct, regParams, badInd] = GetDeformCat3D( expt, infoStruct, varargin) % mouse, exptDate, runs
% Get data from registered 3D calcium imaging experiments
IP = inputParser;
addRequired( IP, 'expt', @isstruct)
addRequired( IP, 'infoStruct', @isstruct)
addOptional( IP, 'deformLim', [], @isstruct)
addParameter( IP, 'zoom', 2, @isnumeric ) 
addParameter( IP, 'edge', nan(1,4), @isnumeric ) %
addParameter( IP, 'basePrct', 10, @isnumeric) %rolling percentile value
addParameter( IP, 'window', 101, @isnumeric) %rolling percentile window
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, infoStruct, varargin{:} ); %expt, 
digiZoom = IP.Results.zoom;
edge = IP.Results.edge;
deformLim = IP.Results.deformLim;
basePrct = IP.Results.basePrct; 
windowSize = IP.Results.window;
if mod(windowSize,2) ==0,  windowSize = windowSize+1; end
show = IP.Results.show;
overwrite = IP.Results.overwrite;
savePath =sprintf('%s%s_deformation.mat', infoStruct.dir, infoStruct.exptName);
regParams = struct();
if ~exist(savePath, 'file') || overwrite
    % Determine limits of real deformation (vs TurboReg failed)
    if isempty(deformLim)
        deformLim = struct('trans',[-Inf, Inf], 'scale',[-Inf, Inf], 'stretch',[-Inf,Inf], 'shear',[-Inf, Inf], 'shift',[-Inf, Inf]); %'stretch',[-Inf,Inf], 
    end
    
    % Initialize structure to store deformation data
    deformResultStruct = struct('RS_final',nan(infoStruct.totScan, infoStruct.Nplane), 'CS_final',nan(infoStruct.totScan, infoStruct.Nplane),'ZS_final',nan(infoStruct.totScan, infoStruct.Nplane),...
        'transAP',nan(infoStruct.totScan, infoStruct.Nplane), 'transML',nan(infoStruct.totScan, infoStruct.Nplane),...
        'scaleAP',nan(infoStruct.totScan, infoStruct.Nplane), 'scaleML',nan(infoStruct.totScan, infoStruct.Nplane),...
        'shearAP',nan(infoStruct.totScan, infoStruct.Nplane), 'shearML',nan(infoStruct.totScan, infoStruct.Nplane),...
        'shiftZ',nan(infoStruct.totScan, infoStruct.Nplane) ); % 'stretchAP',nan(infoStruct.totScan, infoStruct.Nplane), 'stretchML',nan(infoStruct.totScan, infoStruct.Nplane),

    % Get the results of rigid DFT shifting and z interpolation
    if infoStruct.Nplane > 1
        fprintf('\nGetting results of DFT registration');
        %{ 
        %  Old, run-based version
        for runs = flip(1:infoStruct.Nrun)
            [~,shiftPath] = FileFinder( infoStruct.runDir{runs}, 'type','mat', 'contains','dftshifts'); % infoStruct(r).dir
            shiftData(runs) = load(shiftPath{1}, 'RS_final', 'CS_final', 'ZS_final'); %#ok<*AGROW>
        end
        deformResultStruct.RS_final = [shiftData.RS_final]'; 
        deformResultStruct.CS_final = [shiftData.CS_final]';
        deformResultStruct.ZS_final = [shiftData.ZS_final]';
        %}
        [~,shiftPath] = FileFinder( expt.dir, 'type','mat', 'contains','dftshifts');
        if exist(shiftPath{1}, 'file')
            fprintf('\nLoading %s', shiftPath{1})
            shiftData = load(shiftPath{1});
        end
        deformResultStruct.RS_final = [shiftData.RS_final]';
        deformResultStruct.CS_final = [shiftData.CS_final]';
        deformResultStruct.ZS_final = [shiftData.ZS_final]';

        fprintf('\nGetting results of z interpolation');
        % Check if there's a main-folder zinterp file
        [~,interpPath] = FileFinder( expt.dir, 'type','mat', 'contains','zinterp');
        if ~isempty(interpPath)
            interpData = load(interpPath{1});
            deformResultStruct.shiftZ = interpData.ZS'; %[interpData.ZS_chunk]' + [interpData.ZS1]';
        else
            for runs = flip(1:infoStruct.Nrun)
                [~,interpPath] = FileFinder( infoStruct.runDir{runs}, 'type','mat', 'contains','zinterp'); % infoStruct(r).dir
                %interpPath = FileFind( infoStruct.runDir{r}, 'mat', false, @(x)(contains( x, 'zinterp' )) );
                interpData(runs) = load(interpPath{1});
            end
            deformResultStruct.shiftZ = [interpData.ZS_chunk]' + [interpData.ZS1]';
        end
    end

    % Get results of planar affine registration
    fprintf('\nGetting results of planar affine registration');
    regTforms = []; affineParams = [];   
    [~, tformPath] = FileFinder(infoStruct.dir, 'type','mat', 'contains','regTforms'); %FileFind( infoStruct.dir, 'mat', false, @(x)(contains( x, '_affine_tforms' )) );
    if ~isempty(tformPath)
        % Load registration results
        regData = load( tformPath{1} ); %strcat(infoStruct.dir,infoStruct.name,'_affine_tforms.mat')
        if isfield(regData, 'regTform')
            regTforms = regData.regTform;
        elseif isfield(regData, 'tforms_all')
            regTforms = regData.tforms_all;
        end
        % Determine edges used for affine registration
        if isfield(regData, 'params')
            regParams = regData.params;
        end
        if any(isnan(edge))
            if isfield(regData, 'edges')
                edge = regData.edges;
            elseif isfield(regData, 'params') && isfield(regData.params, 'edges')
                edge = regData.params.edges;
            else
                error('Edges unknown');
            end
        end
        fprintf('\nEdges used for affine registration: [L, R, T, B] = [%i, %i, %i, %i]', edge )
        % Get conversion factors
        %umPerPixel = (1/0.53)/digiZoom;
        MLpix = infoStruct(1).sz(1)-edge(3)-edge(4);
        APpix = infoStruct(1).sz(2)-edge(1)-edge(2);
        dT = infoStruct.Nplane/infoStruct.framerate;

        zGood = find(~any( cellfun( @isempty, regTforms ), 2 ))';
        if numel(zGood) < infoStruct.Nplane, warning('affine transform data missing for some planes'); end
        for z = zGood 
            for s = 1:infoStruct.totScan
                deformResultStruct.transAP(s,z) = regTforms{z,s}.T(3,1); % Trans AP (right/left +/-) trans_x
                deformResultStruct.transML(s,z) = regTforms{z,s}.T(3,2); % Trans ML (down/up +/-)  trans_y
                deformResultStruct.scaleAP(s,z) = regTforms{z,s}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
                deformResultStruct.scaleML(s,z) = regTforms{z,s}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
                deformResultStruct.shearAP(s,z) = regTforms{z,s}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
                deformResultStruct.shearML(s,z) = regTforms{z,s}.T(2,1); % Shear ML (tilt down/right +/-) shear_y
            end
        end
        % Calculate stretch here so it can be used to censor bad datapoints
        deformResultStruct.stretchAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(expt.umPerPixel*APpix./deformResultStruct.scaleAP, 1, 1)];
        deformResultStruct.stretchML = (1/dT)*[nan(1,infoStruct.Nplane); diff(expt.umPerPixel*MLpix./deformResultStruct.scaleML, 1, 1)];
    else
        error('\nNo AlignPlanes transform data found');
    end

    % Find and suppress frames with at least one bad deformation value
    badMat = ...
    (deformResultStruct.transAP < deformLim.trans(1) | deformResultStruct.transAP > deformLim.trans(2)) + (deformResultStruct.transML < deformLim.trans(1) | deformResultStruct.transML > deformLim.trans(2)) + ... 
    (deformResultStruct.scaleAP < deformLim.scale(1) | deformResultStruct.scaleAP > deformLim.scale(2)) + (deformResultStruct.scaleML < deformLim.scale(1) | deformResultStruct.scaleML > deformLim.scale(2)) + ... 
    (deformResultStruct.shearAP < deformLim.shear(1) | deformResultStruct.shearAP > deformLim.shear(2)) + (deformResultStruct.shearML < deformLim.shear(1) | deformResultStruct.shearML > deformLim.shear(2)) + ... 
    (deformResultStruct.shiftZ < deformLim.shift(1) | deformResultStruct.shiftZ > deformLim.shift(2));
        %(deformResultStruct.stretchAP < deformLim.stretch(1) | deformResultStruct.stretchAP > deformLim.stretch(2)) + (deformResultStruct.stretchML < deformLim.stretch(1) | deformResultStruct.stretchML > deformLim.stretch(2)) + ... 
    badInd = find(badMat);
    fprintf('\nFound %i bad deformation frames, set to NaN\n', numel(badInd));
    deformCensStruct = deformResultStruct;
    deformCensStruct.transAP(badInd) = NaN;
    deformCensStruct.transML(badInd) = NaN;
    deformCensStruct.scaleAP(badInd) = NaN;
    deformCensStruct.scaleML(badInd) = NaN;
    deformCensStruct.shearAP(badInd) = NaN;
    deformCensStruct.shearML(badInd) = NaN;
    deformCensStruct.shiftZ(badInd) = NaN;

    % Break the results into submovie level structures, and convert to micron units
    deformSubStruct = repmat( struct('dft_RS',[], 'dft_CS',[], 'dft_ZS',[],...
        'transAP',[], 'transML',[], 'transMag',[], 'transAngle',[], 'DtransAP',[], 'DtransML',[], 'DtransMag',[], 'DtransAngle',[],...
        'scaleAP',[], 'scaleML',[], 'scaleMag',[], 'scaleAngle',[], 'stretchAP',[], 'stretchML',[], 'stretchMag',[], 'stretchAngle',[],...
        'shearAP',[], 'shearML',[], 'shearMag',[], 'shearAngle',[], 'DshearAP',[], 'DshearML',[], 'DshearMag',[], 'DshearAngle',[],... 
        'shiftZ',[], 'compressZ',[], 'DshiftZ',[], 'DcompressZ',[]), 1, infoStruct.Nrun );
    scanLims = [0, cumsum(infoStruct.Nscan)];

    for runs = 1:infoStruct.Nrun
        subScans = scanLims(runs)+1:scanLims(runs+1);
        % DFT REGISTRATION RESULTS
        deformSubStruct(runs).dft_RS = expt.umPerPixel*deformCensStruct.RS_final(subScans,:);
        deformSubStruct(runs).dft_CS = expt.umPerPixel*deformCensStruct.CS_final(subScans,:);
        deformSubStruct(runs).dft_ZS = deformCensStruct.ZS_final(subScans);
        
        % Z INTERPOLATION RESULTS
        deformSubStruct(runs).shiftZ = deformCensStruct.shiftZ(subScans,:);
        
        % AFFINE REGISTRATION RESULTS (converting to um in the process)
        deformSubStruct(runs).transAP = expt.umPerPixel*deformCensStruct.transAP(subScans,:);
        deformSubStruct(runs).transML = expt.umPerPixel*deformCensStruct.transML(subScans,:);
        deformSubStruct(runs).scaleAP = expt.umPerPixel*APpix./deformCensStruct.scaleAP(subScans,:); %expt.umPerPixel*APpix*deformCatStruct.scaleAP(subScans,:);
        deformSubStruct(runs).scaleML = expt.umPerPixel*MLpix./deformCensStruct.scaleML(subScans,:); %expt.umPerPixel*MLpix*deformCatStruct.scaleML(subScans,:);
        deformSubStruct(runs).shearAP = deformCensStruct.shearAP(subScans,:);
        deformSubStruct(runs).shearML = deformCensStruct.shearML(subScans,:);

        % High-pass filter (optional)
        if basePrct > 0
            fprintf('\nHP Filtering: window = %i scans, %ith percentile subtraction', windowSize, basePrct);
            deformSubStruct(runs).transAP = deformSubStruct(runs).transAP - MovingPercentile(deformSubStruct(runs).transAP, basePrct, windowSize, 'pre');
            deformSubStruct(runs).transML = deformSubStruct(runs).transML - MovingPercentile(deformSubStruct(runs).transML, basePrct, windowSize, 'pre');
            deformSubStruct(runs).scaleAP = deformSubStruct(runs).scaleAP - MovingPercentile(deformSubStruct(runs).scaleAP, basePrct, windowSize, 'pre');
            deformSubStruct(runs).scaleML = deformSubStruct(runs).scaleML - MovingPercentile(deformSubStruct(runs).scaleML, basePrct, windowSize, 'pre');
            deformSubStruct(runs).shearAP = deformSubStruct(runs).shearAP - MovingPercentile(deformSubStruct(runs).shearAP, basePrct, windowSize, 'pre');
            deformSubStruct(runs).shearML = deformSubStruct(runs).shearML - MovingPercentile(deformSubStruct(runs).shearML, basePrct, windowSize, 'pre');
            deformSubStruct(runs).shiftZ = deformSubStruct(runs).shiftZ - MovingPercentile(deformSubStruct(runs).shiftZ, basePrct, windowSize, 'pre');
        end
        % Calculate derivatives
        deformSubStruct(runs).compressZ = [nan(numel(subScans),1), diff( deformCensStruct.shiftZ(subScans,:), 1, 2 )]; % spatial derivative
        deformSubStruct(runs).DtransAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(runs).transAP, 1, 1)];
        deformSubStruct(runs).DtransML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(runs).transML, 1, 1)];
        deformSubStruct(runs).stretchAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(runs).scaleAP, 1, 1)];
        deformSubStruct(runs).stretchML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(runs).scaleML, 1, 1)];
        deformSubStruct(runs).DshearAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(runs).shearAP, 1, 1)];
        deformSubStruct(runs).DshearML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(runs).shearML, 1, 1)];
        deformSubStruct(runs).DshiftZ = (1/dT)*[nan(1,infoStruct.Nplane); diff( deformSubStruct(runs).shiftZ, 1, 1)];
        deformSubStruct(runs).DcompressZ = (1/dT)*[nan(1,infoStruct.Nplane); diff( deformSubStruct(runs).compressZ, 1, 1)];
        % Convert xy variables to polar coordinates
        [deformSubStruct(runs).transAngle, deformSubStruct(runs).transMag] = cart2pol( deformSubStruct(runs).transAP, deformSubStruct(runs).transML ); 
        [deformSubStruct(runs).scaleAngle, deformSubStruct(runs).scaleMag] = cart2pol( deformSubStruct(runs).scaleAP, deformSubStruct(runs).scaleML ); 
        [deformSubStruct(runs).shearAngle, deformSubStruct(runs).shearMag] = cart2pol( deformSubStruct(runs).shearAP, deformSubStruct(runs).shearML ); 
        [deformSubStruct(runs).DtransAngle, deformSubStruct(runs).DtransMag] = cart2pol( deformSubStruct(runs).DtransAP, deformSubStruct(runs).DtransML );
        [deformSubStruct(runs).stretchAngle, deformSubStruct(runs).stretchMag] = cart2pol( deformSubStruct(runs).stretchAP, deformSubStruct(runs).stretchML );
        [deformSubStruct(runs).DshearAngle, deformSubStruct(runs).DshearMag] = cart2pol( deformSubStruct(runs).DshearAP, deformSubStruct(runs).DshearML);
        % Break up expansive, compressive, and mixed stretch
        expInd = deformSubStruct(runs).stretchAngle >= 0 & deformSubStruct(runs).stretchAngle <= pi/2;
        compInd = deformSubStruct(runs).stretchAngle <= -pi/2 & deformSubStruct(runs).stretchAngle >= -pi;
        mixedInd = (deformSubStruct(runs).stretchAngle < pi & deformSubStruct(runs).stretchAngle > pi/2) | (deformSubStruct(runs).stretchAngle < 0 & deformSubStruct(runs).stretchAngle > -pi/2);
        deformSubStruct(runs).stretchExp = zeros( size(deformSubStruct(runs).scaleAP) );
        deformSubStruct(runs).stretchExp( expInd ) = deformSubStruct(runs).stretchMag( expInd );
        deformSubStruct(runs).stretchComp = zeros( size(deformSubStruct(runs).scaleAP) );
        deformSubStruct(runs).stretchComp( compInd ) = deformSubStruct(runs).stretchMag( compInd );
        deformSubStruct(runs).stretchMixed = zeros( size(deformSubStruct(runs).scaleAP) );
        deformSubStruct(runs).stretchMixed( mixedInd ) = deformSubStruct(runs).stretchMag( mixedInd );

    end
 
    fprintf('\nSaving %s', savePath);
    save(savePath, 'infoStruct','deformLim','deformSubStruct','deformResultStruct','deformCensStruct','badInd','expInd','compInd','mixedInd',...
        'affineParams','tformPath','dT','APpix','MLpix','edge','basePrct','windowSize','digiZoom','zGood','scanLims'); % ,'umPerPixel'
else
    fprintf('\nLoading %s', savePath);
    load(savePath);
end

if show
    runTicks = cumsum(infoStruct.Nscan);
    timeTicks = 1:5000:infoStruct.totScan;
    TL = [0.005,0];
    leftOpt = {[0.01,0.09], [0.06,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right]}
    rightOpt = {[0.01,0.1], [0.06,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right]}
    close all; clearvars sp
    figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
    if infoStruct.Nplane == 1
        planeTicks = 1;
        % Translation 
        sp(1) = subtightplot(6, 2, 1, rightOpt{:}); 
        MakeDeformPlot( deformResultStruct.transAP, 'AP Translation', TL, runTicks, planeTicks ) % , deformLim.trans
        hold on;
        line([1,infoStruct.totScan], deformLim.trans(1)*[1,1], 'color','r', 'linestyle','--')
        line([1,infoStruct.totScan], deformLim.trans(2)*[1,1], 'color','r', 'linestyle','--')
        title('Turboreg Registration Results');
        sp(2) = subtightplot(6, 2, 3, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transML, 'ML Translation', TL, runTicks, planeTicks ) % , deformLim.trans
        hold on;
        line([1,infoStruct.totScan], deformLim.trans(1)*[1,1], 'color','r', 'linestyle','--')
        line([1,infoStruct.totScan], deformLim.trans(2)*[1,1], 'color','r', 'linestyle','--')
        % Scaling
        sp(3) = subtightplot(6, 2, 5, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleAP, 'AP Scale', TL, runTicks, planeTicks ) % , deformLim.scale
        hold on;
        line([1,infoStruct.totScan], deformLim.scale(1)*[1,1], 'color','r', 'linestyle','--')
        line([1,infoStruct.totScan], deformLim.scale(2)*[1,1], 'color','r', 'linestyle','--')
        sp(4) = subtightplot(6, 2, 7, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleML, 'ML Scale', TL, runTicks, planeTicks ) % , deformLim.scale
        hold on;
        line([1,infoStruct.totScan], deformLim.scale(1)*[1,1], 'color','r', 'linestyle','--')
        line([1,infoStruct.totScan], deformLim.scale(2)*[1,1], 'color','r', 'linestyle','--')
        % Stretch
        %{
        sp(5) = subtightplot(6, 2, 9, rightOpt{:});
        MakeDeformPlot( deformResultStruct.stretchAP, 'AP Stretch', TL, runTicks, planeTicks, deformLim.stretch )
        sp(6) = subtightplot(6, 2, 11, rightOpt{:});
        MakeDeformPlot( deformResultStruct.stretchML, 'ML Stretch', TL, timeTicks, planeTicks, deformLim.stretch )
        %}
        % Shearing
        sp(5) = subtightplot(6, 2, 9, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearAP, 'AP Shear', TL, runTicks, planeTicks ) % , deformLim.shear
        hold on;
        line([1,infoStruct.totScan], deformLim.shear(1)*[1,1], 'color','r', 'linestyle','--')
        line([1,infoStruct.totScan], deformLim.shear(2)*[1,1], 'color','r', 'linestyle','--')
        sp(6) = subtightplot(6, 2, 11, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearML, 'ML Shear', TL, timeTicks, planeTicks ) % runTicks , deformLim.shear
        hold on;
        line([1,infoStruct.totScan], deformLim.shear(1)*[1,1], 'color','r', 'linestyle','--')
        line([1,infoStruct.totScan], deformLim.shear(2)*[1,1], 'color','r', 'linestyle','--')
        xlabel('Frame');
        
        % Censored final results
        % Translation
        sp(9) = subtightplot(6, 2, 2, rightOpt{:});
        MakeDeformPlot(vertcat(deformSubStruct.transAP), 'AP Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.trans
        title('Final, Censored Deformation');
        sp(10) = subtightplot(6, 2, 4, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transML), 'ML Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Scaling
        sp(11) = subtightplot(6, 2, 6, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleAP), 'AP Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.scale
        sp(12) = subtightplot(6, 2, 8, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleML), 'ML Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] )  
        % Stretch
        %{
        sp(13) = subtightplot(6, 2, 10, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchAP), 'AP Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.stretch
        sp(14) = subtightplot(6, 2, 12, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchML), 'ML Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] )
        %}
        % Shearing
        sp(13) = subtightplot(6, 2, 10, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearAP), 'AP Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        sp(14) = subtightplot(6, 2, 12, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearML), 'ML Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        
        linkaxes(sp,'x');
        xlabel('Frame');
        xlim([1,infoStruct.totScan]) %axis tight;
    else
        if infoStruct.Nplane == 15
            planeTicks = [1:3:15, 15];
        elseif infoStruct.Nplane == 30
            planeTicks = [1:5:30, 30];
        else
            planeTicks = round(linspace(1, infoStruct.Nplane, infoStruct.Nplane/4));
        end
        % Results of preliminary registration
        SP(3) = subtightplot(3,3,7,leftOpt{:});
        MakeDeformPlot( deformResultStruct.ZS_final, 'DFT Z shift', TL, runTicks, planeTicks, deformLim.shift )
        % Row
        SP(1) = subtightplot(3,3,1,leftOpt{:});
        MakeDeformPlot( deformResultStruct.RS_final, 'DFT Row shift', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.shift
        title('Rigid DFT Pre-registration');
        % Column
        SP(2) = subtightplot(3,3,4,leftOpt{:});
        MakeDeformPlot( deformResultStruct.CS_final, 'DFT Column shift', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.shift
        linkaxes(SP,'x');
        axis tight;
        
        % Uncensored Results of affine registration
        % Z Shift
        sp(9) = subtightplot(7, 3, 20, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shiftZ, 'Z Shift', TL, runTicks, planeTicks, deformLim.shift )
        xlabel('Volume Scan');
        % Translation 
        sp(1) = subtightplot(7, 3, 2, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transAP, 'AP Translation', TL, runTicks, planeTicks, deformLim.trans )
        title('Turboreg Registration Results');
        sp(2) = subtightplot(7, 3, 5, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transML, 'ML Translation', TL, runTicks, planeTicks, deformLim.trans )
        % Scaling
        sp(3) = subtightplot(7, 3, 8, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleAP, 'AP Scale', TL, runTicks, planeTicks, deformLim.scale )
        sp(4) = subtightplot(7, 3, 11, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleML, 'ML Scale', TL, runTicks, planeTicks, deformLim.scale ) 
%         % Stretch
%         sp(5) = subtightplot(7, 3, 9, rightOpt{:});
%         MakeDeformPlot( deformResultStruct.stretchAP, 'AP Stretch', TL, runTicks, planeTicks, deformLim.stretch )
%         sp(6) = subtightplot(7, 3, 12, rightOpt{:});
%         MakeDeformPlot( deformResultStruct.stretchML, 'ML Stretch', TL, timeTicks, planeTicks, deformLim.stretch )
%         % Shearing
        sp(7) = subtightplot(7, 3, 14, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearAP, 'AP Shear', TL, runTicks, planeTicks, deformLim.shear )
        sp(8) = subtightplot(7, 3, 17, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearML, 'ML Shear', TL, runTicks, planeTicks, deformLim.shear )

        % Censored final results
        % Z shift
        sp(18) = subtightplot(9, 3, 27, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shiftZ), 'Z Shift (planes)', TL, runTicks, planeTicks, deformLim.shift )
        xlabel('Volume Scan');
        % Translation
        sp(10) = subtightplot(9, 3, 3, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transAP), 'AP Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.trans
        title('Final, Censored Deformation');
        sp(11) = subtightplot(9, 3, 6, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transML), 'ML Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Scaling
        sp(12) = subtightplot(9, 3, 9, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleAP), 'AP Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.scale
        sp(13) = subtightplot(9, 3, 12, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleML), 'ML Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] )  
        % Stretch
        sp(14) = subtightplot(9, 3, 15, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchAP), 'AP Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.stretch
        sp(15) = subtightplot(9, 3, 18, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchML), 'ML Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Shearing
        sp(16) = subtightplot(9, 3, 21, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearAP), 'AP Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        sp(17) = subtightplot(9, 3, 24, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearML), 'ML Shear', TL, runTicks, planeTicks, [-Inf,Inf] )

        linkaxes(sp,'xy');
        xlim([1,infoStruct.totScan]) %axis tight;
        impixelinfo;
        end
end
end

    %{
    % Zposition, as measured by knobby, for Pollen 221130 FOV1
    z_knob = nan(expt.totScan, 1);
    z_knob(expt.scanLims(1)+1:expt.scanLims(2)) = 0; % run 1
    z_knob(expt.scanLims(2)+1:expt.scanLims(3)) = -6.25; % run 2
    z_knob(expt.scanLims(3)+1:expt.scanLims(3)+30) = -6.25; % run 3, phase 1
    z_knob(expt.scanLims(3)+31:expt.scanLims(4)) = -3.9; % run 3, phase 2
    z_knob(expt.scanLims(4)+1:expt.scanLims(4)+30) = -3.9; % run 3, phase 1
    z_knob(expt.scanLims(4)+31:expt.scanLims(5)) = 12.8; % 

    figure;
    subplot(4,1,1); line([1,expt.totScan],[0,0],'color','k'); hold on;
    plot(deformResultStruct.ZS_final); ylabel('ZS-final'); hold on;
    set(gca,'Xtick',expt.scanLims(1:expt.Nruns)+1)
    subplot(4,1,2); line([1,expt.totScan],[0,0],'color','k'); hold on;
    plot(-mean(deformResultStruct.shiftZ,1)'); ylabel('shiftZ'); hold on;
    line([1,expt.totScan],[0,0],'color','k')
    set(gca,'Xtick',expt.scanLims(1:expt.Nruns)+1)
    subplot(4,1,3); line([1,expt.totScan],[0,0],'color','k'); hold on;
    plot(deformResultStruct.ZS_final+ mean(deformResultStruct.shiftZ,1)'); ylabel('sum'); hold on;
    set(gca,'Xtick',expt.scanLims(1:expt.Nruns)+1)
    subplot(4,1,4);  line([1,expt.totScan],[0,0],'color','k'); hold on;
    plot(z_knob); ylabel('z (um)'); hold on;
    set(gca,'Xtick',expt.scanLims(1:expt.Nruns)+1)
    %}