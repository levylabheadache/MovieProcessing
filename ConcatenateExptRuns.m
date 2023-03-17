function interRunShift = ConcatenateExptRuns(expt, runInfo, catInfo, varargin )
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'runInfo', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
addParameter( IP, 'minInt', 1000, @isnumeric )
addParameter( IP, 'maxEdge', [100, 100, 120, 120], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'setEdge', [], @isnumeric )
addParameter( IP, 'refRun', 1, @isnumeric )
addParameter( IP, 'refChan', 'green', @ischar )
addParameter( IP, 'sbx', 'sbxz', @ischar )
addParameter( IP, 'chunkSize', 100, @isnumeric )
addParameter( IP, 'ext', 'sbxcat', @ischar )
addParameter( IP, 'overwrite', false, @islogical ) % for scanbox, 1 = green, 2 = red. -1 = both
parse( IP, expt, runInfo, catInfo, varargin{:} ); % mouse, exptDate,
minIntOrig = IP.Results.minInt;
if expt.Nchan == 2 && expt.Nplane > 1
    minInt = minIntOrig/2^8;
else
    minInt = minIntOrig;
end
edgeMax = IP.Results.maxEdge;
edgeSet = IP.Results.setEdge;
sbxType = IP.Results.sbx;
refRun = IP.Results.refRun;
refChan = IP.Results.refChan;
refChanInd = find(strcmpi(refChan, {'red','green'} ));
if isempty(refChanInd),  error('Invalid reference channel'); end
catExt = IP.Results.ext; %'.sbxcat'; % '.sbxcat'
chunkSize = IP.Results.chunkSize;
overwrite = IP.Results.overwrite;

catName = expt.name; % sprintf('%s_FOV%i', , expt.fov);
catPathRoot = sprintf('%s%s', expt.dir, catName);
catSbxPath = strcat(catPathRoot, '.', catExt);
interRunShift = zeros(expt.Nruns, 3);
if expt.Nruns > 1
    catDir = strcat(expt.dir,'Concat\'); mkdir(catDir)
    if overwrite || ~exist(catSbxPath,'file')
        Nchan = runInfo(1).nchan; Nscan = [runInfo.Nscan];  %totScan = sum(Nscan);
        sbxPath = cell(1,expt.Nruns);  runProjPath = cell(1,expt.Nruns);  runProj = cell(expt.Nruns, 1);
        cropProj = cell(expt.Nruns, 1); shiftProj = cell(1,expt.Nruns); nullProj = cell(1,expt.Nruns); %hpProj = cell(expt.Nruns, 1);
        if isempty(edgeSet), projEdges = zeros( expt.Nruns, 4 ); else, projEdges = repmat(edgeSet, expt.Nruns, 1 ); end
        if expt.Nplane > 1
            Nplane = runInfo(1).Nplane;   % Nrow = runInfo(1).sz(1); Ncol = runInfo(1).sz(2);
            for runs = 1:expt.Nruns % expt.runs
                sbxPath{runs} = sprintf('%s%s.%s', runInfo(runs).dir, runInfo(runs).fileName, sbxType ); % sbx_affine  sbxz
                runProj{runs} = WriteSbxProjection(sbxPath{runs}, runInfo(runs), 'verbose',true, 'chan',refChan, 'monochrome',true, 'RGB',false, 'type','dft', 'overwrite',false); % 'Z'
                if ndims(runProj{runs}) == 4, runProj{runs} = squeeze(runProj{runs}(:,:,refChanInd,:)); end
                if isempty(edgeSet)
                    projEdges(runs,:) = GetEdges( runProj{runs}(:,:,round(expt.Nplane/2)), 'minInt',minInt, 'show',true );
                else
                    ShowEdges(projEdges(runs,:), runProj{runs}(:,:,round(expt.Nplane/2)));
                    %pause;
                end
                %saveastiff( runProj{r}, sprintf('%sTest\\%s_Proj_run%i.tif', expt.dir, expt.name, r ) );
            end
            close all;
            % keep all edges within edgeMax
            if isempty(edgeSet)
                edgeSet = max(projEdges, [], 1); %[90,100,60,100]; % %edgeSet(2) =  edgeSet(2) + (683-539)
            end
            aboveEdge = edgeMax - edgeSet < 0;
            edgeSet(aboveEdge) = edgeMax(aboveEdge);

            %filtSigma = 5;
            for runs = 1:expt.Nruns % expt.runs
                cropProj{runs} = runProj{runs}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:);
                %hpProj{r} = cropProj{r} - imgaussfilt(cropProj{r}, filtSigma);
                %saveastiff( uint16(cropProj{runs}), sprintf('%s%s_Proj_run%i.tif', catDir, expt.name, runs ) );
                %saveastiff( hpProj{r}, sprintf('%s%s_HP_run%i.tif', catDir, expt.name, r ) );
            end
            % Save the uncorrected projections as a BioFormats tif
            cropProjCat = cat( 4, cropProj{:} );
            bfsave(cropProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:,:), sprintf('%s%s_uncorrected.tif', catDir, expt.name ))

            % Calculate shifts between runs
            zUse = 1:Nplane;
            %refRun = 2;  %10:Nplane;
            interRunDFTshift = zeros(expt.Nruns,3);  interRunZshift = zeros(numel(zUse), expt.Nruns);
            for runs = 1:expt.Nruns
                if runs ~= refRun
                    if Nplane > 1
                        interRunDFTshift(runs,:) = dftregistration3D(fftn(cropProj{refRun}(:,:,zUse)),  fftn(cropProj{runs}(:,:,zUse)), 4);
                        interRunZshift(:,runs) = InterpolateZshift(cropProj{refRun}(:,:,zUse),  cropProj{runs}(:,:,zUse), 2);
                    else
                        zUse = 1;
                        tempOut = dftregistration(fft2(cropProj{refRun}),  fft2(cropProj{runs})); % [error,diffphase,net_row_shift,net_col_shift]
                        interRunDFTshift(runs,1:2) = tempOut(3:4);
                    end
                end
            end
            %interRunShift = zeros(expt.Nruns,3);
            interRunShift(:,[1,2,3]) = interRunDFTshift(:,[2,1,3]); % imtranslate expects [x,y,z], not [y,x,z]
            interRunShift = round( interRunShift );

            % Apply shifts to each run's projection
            shiftProj = cell(expt.Nruns,expt.Nchan);
            for runs = 1:expt.Nruns
                if Nplane > 1
                    shiftProj{runs} = imtranslate( runProj{runs}, interRunShift(runs,:) );
                else
                    shiftProj{runs} = imtranslate( runProj{runs}, interRunShift(runs,1:2) );
                end
                %saveastiff( [cropProj{refRun}, shiftProj{runs}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:)], sprintf('%s%s_shiftProj_run%i.tif', catDir, expt.name, runs ) );
            end
            shiftProjCat = cat( 4, shiftProj{:} ); % projCat = cat( 4, nullProj{:} );
            bfsave(shiftProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:,:), sprintf('%s%s_cat.tif', catDir, expt.name ))
        else
            for runs = expt.runs
                sbxPath{runs} = runInfo(runs).path; %sprintf('%s%s.sbx', runInfo(r).dir, runInfo(r).fileName );
                runProj{runs} = WriteSbxProjection(runInfo(runs).path, runInfo(runs), 'binT',1, 'verbose',true, 'chan',refChan, 'monochrome',true, 'type','dft', 'overwrite',false);
                if isempty(edgeSet)
                    projEdges(runs,:) = GetEdges( runProj{runs}(:,:,end), 'minInt',minInt, 'show',true );
                else
                    ShowEdges(projEdges(runs,:), runProj{runs});
                end
            end
            close all;

            % keep all edges within edgeMax
            if isempty(edgeSet)
                edgeSet = max(projEdges, [], 1); %[90,100,60,100]; % %edgeSet(2) =  edgeSet(2) + (683-539)
            end
            aboveEdge = edgeMax - edgeSet < 0;
            edgeSet(aboveEdge) = edgeMax(aboveEdge);
            % Crop the reference projections
            for runs = expt.runs
                cropProj{runs} = runProj{runs}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:);
                %saveastiff( cropProj{runs}, sprintf('%s%s_Proj_run%i.tif', catDir, expt.name, runs ) );
            end
            close all;

            % Calculate shifts between runs
            refFT = fft2(cropProj{refRun});
            interRunDFTshift = zeros(expt.Nruns,2);
            for runs = 1:expt.Nruns
                if runs ~= refRun
                    tempOut = dftregistration(refFT,  fft2(cropProj{runs})); % [error,diffphase,net_row_shift,net_col_shift]  hpProj
                    interRunDFTshift(runs,1:2) = tempOut(3:4);
                end
            end
            interRunShift = round( flip(interRunDFTshift, 2) ); % flip shifts to be compatible with imtranslate

            % Apply shifts to each run's projection
            for runs = 1:expt.Nruns
                nullProj{runs} = imtranslate( runProj{runs}, interRunShift(runs,1:2) );
                %saveastiff( nullProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:), sprintf('%s%s_nullProj_run%i.tif', catDir, expt.name, r ) );
                shiftProj{runs} = nullProj{runs};
                %saveastiff( [cropProj{refRun}, shiftProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:)], sprintf('%s%s_shiftProj_run%i.tif', catDir, expt.name, r ) );
            end
            projCat = uint16(cat( 4, nullProj{:} )) ;
            shiftProjCat = uint16(cat( 4, shiftProj{:} ));
            for z = 1:expt.Nplane
                saveastiff( squeeze( projCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_ind_z%i.tif', catDir, expt.name, z ) )
                saveastiff( squeeze( shiftProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_cat_z%i.tif', catDir, expt.name, z ) )
            end
        end

        % Write the sbxcat file
        fprintf('\n     Writing %s\n', catSbxPath); tic
        rw = SbxWriter(catSbxPath, catInfo, catExt, true); 
        w = waitbar(0, sprintf('writing %s',catExt));
        for runs = 1:expt.Nruns %expt.runs
            % Load the run stack, in chunks
            [runChunkLims, NrunChunk, runChunkLength] = MakeChunkLims(1, runInfo(runs).Nscan, runInfo(runs).Nscan, 'allowPartial',true, 'size',chunkSize);
            fprintf('\n   Loading %s (%i chunks at a time)... ', sbxPath{runs}, chunkSize); tic
            if expt.Nplane > 1
                for chunk = 1:NrunChunk
                    runChunkStack = readSBX(sbxPath{runs}, runInfo(runs), runChunkLims(chunk,1), runChunkLength(chunk), -1, []); % [c,x,y,z,t]
                    if Nchan == 2
                        runChunkStack = permute(runChunkStack, [2,3,4,5,1]); % [x,y,z,t,c]
                        % Apply shifts to each channel of each scan
                        if ~all(interRunShift(runs,:) == 0)
                            for s = 1:size(runChunkStack, 4) %1:Nscan(runs)
                                for c = 1:Nchan
                                    runChunkStack(:,:,:,s,c) = imtranslate( runChunkStack(:,:,:,s,c), interRunShift(runs,:)  );
                                end
                            end
                        end
                        runChunkStack = permute(runChunkStack, [5,1,2,3,4]); % [c,x,y,z,t]
                        runChunkStack = reshape(runChunkStack, [size(runChunkStack,[1,2,3]), prod(size(runChunkStack,[4,5]))]); % rw expects this form
                    elseif Nchan == 1
                        % Apply shifts to each scan
                        if ~all(interRunShift(runs,:) == 0)
                            for s = 1:size(runChunkStack, 4) 
                                runChunkStack(:,:,:,s) = imtranslate( runChunkStack(:,:,:,s), interRunShift(runs,:)  );
                            end
                        end
                        runChunkStack = reshape(runChunkStack, [size(runChunkStack,[1,2]), prod(size(runChunkStack,[3,4]))] );
                    end
                    rw.write( runChunkStack ); % Write the chunk to sbxcat
                end
            else
                for chunk = 1:NrunChunk
                    runChunkStack = readSBX(sbxPath{runs}, runInfo(runs), runChunkLims(chunk,1), runChunkLength(chunk), -1, []); % [c,x,y,t]
                    if Nchan == 2
                        runChunkStack = permute(runChunkStack, [2,3,4,1]); % [x,y,t,c]
                        % Apply shifts to each color of each scan
                        if ~all(interRunShift(runs,:) == 0)
                            for s = 1:size(runChunkStack, 3) 
                                for c = 1:Nchan
                                    runChunkStack(:,:,s,c) = imtranslate( runChunkStack(:,:,s,c), interRunShift(runs,:)  );
                                end
                            end
                        end
                        runChunkStack = permute(runChunkStack, [4,1,2,3]); % [c,x,y,t]
                    elseif Nchan == 1
                        % Apply shifts to each scan
                        if ~all(interRunShift(runs,:) == 0)
                            for s = 1:size(runChunkStack, 3) 
                                runChunkStack(:,:,s) = imtranslate( runChunkStack(:,:,s), interRunShift(runs,:)  );
                            end
                        end
                    end
                    rw.write( runChunkStack ); % Write the chunk to sbxcat
                end
            end
            waitbar( runs/expt.Nruns, w );
            toc
        end
        rw.delete;
        delete(w);
    else
        fprintf('\n%s already exists!', catSbxPath);
    end
elseif ~exist(catSbxPath, 'file') || overwrite
    fprintf('\nSingle run experiment: copying %s to %s', runInfo.path, catSbxPath);
    copyfile(runInfo.path, catSbxPath);
    %[~,~,projPath] = WriteSbxProjection(runInfo.path, runInfo, 'monochrome',true', 'RGB',true );
    %{
    % copy-paste the projections too
    if expt.Nplane == 1
        FileFinder(runInfo.dir, 'contains','')
    else
        FileFinder(runInfo.dir, 'contains','_Z_')
    end
    %}
    %projSourcePath = strcat(runInfo.dir, runInfo.fileName, '_interpProj.tif');
    %copyfile(projSourcePath, catProjPath);
end
end