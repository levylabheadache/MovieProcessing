function interRunShift = ConcatenateSinglePlane(expt, runInfo, catInfo, refRun, commonEdge)

catExt = '.sbxcat'; %'.sbxz'; %
sbxCatPath = [expt.dir, expt.name, catExt];   % [expt.dir, expt.name,'_redo', catExt];
if ~exist(sbxCatPath,'file')
    if expt.Nruns > 1
        sbxPath = cell(1,expt.Nruns);  projPath = cell(1,expt.Nruns);
        runProj = cell(expt.Nruns, 1); cropProj = cell(expt.Nruns, 1); hpProj = cell(expt.Nruns, 1);
        mkdir(strcat(expt.dir,'Test'))
        for r = expt.runs
            sbxPath{r} = sprintf('%s%s.sbx', runInfo(r).dir, runInfo(r).fileName ); % sbx_affine  sbxz
            projPath{r} = sprintf('%s%s_rawProj_green.tif', runInfo(r).dir, runInfo(r).fileName ); % affine
            runProj{r} = loadtiff( projPath{r} );
        end
        close all;

        filtSigma = 5;
        for r = expt.runs
            cropProj{r} = runProj{r}(commonEdge(3):end-commonEdge(4),commonEdge(1):end-commonEdge(2),:);
            hpProj{r} = cropProj{r} - imgaussfilt(cropProj{r}, filtSigma);
            saveastiff( cropProj{r}, sprintf('%sTest\\%s_Proj_run%i.tif', expt.dir, expt.name, r ) );
            saveastiff( hpProj{r}, sprintf('%sTest\\%s_HP_run%i.tif', expt.dir, expt.name, r ) );
        end

        % Calculate shifts between runs
        interRunDFTshift = zeros(expt.Nruns,2);
        for r = 1:expt.Nruns
            if r ~= refRun
                tempOut = dftregistration(fft2(cropProj{refRun}),  fft2(cropProj{r})); % [error,diffphase,net_row_shift,net_col_shift]  hpProj
                interRunDFTshift(r,1:2) = tempOut(3:4);
            end
        end
        interRunShift = round( flip(interRunDFTshift, 2) ); % flip shifts to be compatible with imtranslate

        % Apply shifts to each run's projection
        shiftProj = cell(1,expt.Nruns); nullProj = cell(1,expt.Nruns);
        for r = 1:expt.Nruns
            nullProj{r} = imtranslate( runProj{r}, interRunShift(r,1:2) );
            saveastiff( nullProj{r}(commonEdge(3):end-commonEdge(4),commonEdge(1):end-commonEdge(2),:), sprintf('%sTest\\%s_nullProj_run%i.tif', expt.dir, expt.name, r ) );
            if expt.Nplane > 1
                shiftProj{r} = imtranslate( runProj{r}, interRunShift(r,:) );
            else
                shiftProj{r} = nullProj{r};
            end
            saveastiff( [cropProj{refRun}, shiftProj{r}(commonEdge(3):end-commonEdge(4),commonEdge(1):end-commonEdge(2),:)], sprintf('%sTest\\%s_shiftProj_run%i.tif', expt.dir, expt.name, r ) );
        end
        projCat = cat( 4, nullProj{:} ); shiftProjCat = cat( 4, shiftProj{:} );
        for z = 1:expt.Nplane
            saveastiff( squeeze( projCat(commonEdge(3):end-commonEdge(4),commonEdge(1):end-commonEdge(2),z,:) ), sprintf('%sTest\\%s_ind_z%i.tif', expt.dir, expt.name, z ) )
            saveastiff( squeeze( shiftProjCat(commonEdge(3):end-commonEdge(4),commonEdge(1):end-commonEdge(2),z,:) ), sprintf('%sTest\\%s_cat_z%i.tif', expt.dir, expt.name, z ) )
        end

        % Write the concatenated sbx file
        fprintf('\n     Writing %s\n', sbxCatPath); tic
        w = waitbar(0, sprintf('writing %s',catExt));
        rw = SbxWriter(sbxCatPath, catInfo, catExt, true); %pipe.io.RegWriter(sbxCatPath, catInfo, catExt, true); % info,
        for r = expt.runs % %
            % Load run stack
            fprintf('\n   Loading %s... ', sbxPath{r}); tic
            runStack = readSBX(sbxPath{r}, runInfo(r), 1, runInfo(r).Nscan, -1, []); toc
            runStack = permute( runStack, [2,3,4,1]);
            runStack = imtranslate(runStack, interRunShift(r,1:2));
            runStack = permute( runStack, [4,1,2,3]);
            rw.write( runStack ); %
            waitbar( r/expt.Nruns, w );
            toc
        end
        rw.delete;
        delete(w);
    elseif expt.Nruns == 1
        copyfile(runInfo.path, sbxCatPath)
    end
else
    interRunShift = nan(expt.Nruns, 2);
    fprintf('\n%s already exists  ', sbxCatPath);
end

% Write concatenated, binned tif and mean projection
catProjPath = strcat(expt.dir, expt.name,'_catProj.tif');
if ~exist(catProjPath, 'file') %true %
    catStack = WriteSbxPlaneTif(sbxCatPath, catInfo, 1, 'dir',strcat(catInfo.dir,'Ztifs\'), 'name',expt.name, 'type','raw', 'edge',[0,0,0,0], ...
        'scale',1, 'binT',10, 'verbose',true, 'overwrite',false ); %, 'chan',refChan   , 'Nscan',55800 , 'Nscan',100, 'firstScan',55801
    catMean = mean( catStack, 3);
    saveastiff(uint16(catMean), catProjPath );
end
end