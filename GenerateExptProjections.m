function projParam = GenerateExptProjections(expt, catInfo, Tscan, varargin) 

% Generate downsampled, possibly z-projected, movies for each run from each requested type of concatenated data
if isempty(varargin)
    projParam.overwrite = false;
else
    projParam = varargin{1};
end
projParam.dir = strcat(expt.dir, 'Projections\'); 
mkdir(projParam.dir);
projParamPath = sprintf('%s%s_projParam.mat', projParam.dir, expt.name);
chanName = {'red', 'green'}; 
if nargin > 3  %~exist(projParamPath, 'file') % || projParam.overwrite
    runsDir = sprintf('%sRuns\\',projParam.dir); mkdir(runsDir);
    if expt.Nplane == 1
       projParam.z = {1}; 
    end
    projParam.Ncolor = numel(projParam.color);
    projParam.Nz = numel(projParam.z);
    projParam.bin = max([1, round(expt.scanRate/projParam.rate_target)]);
    projParam.rate_bin = expt.scanRate/projParam.bin;
    projParam.scaleFactor = projParam.umPerPixel_target/expt.umPerPixel; % round(projParam.umPerPixel_target/expt.umPerPixel); %
    if projParam.scaleFactor < 1
       projParam.scaleFactor = 1; 
    end
    projParam.umPerPixel_scale = expt.umPerPixel*projParam.scaleFactor;
    disp(projParam);

    fprintf('\nGenerating projections:');
    Tcat = vertcat(Tscan{:});
    projParam.Ntype = numel(projParam.sbx_type);
    for t = 1:projParam.Ntype
        if exist(expt.sbx.(projParam.sbx_type{t}), 'file')
            fprintf('\n  %s', projParam.sbx_type{t});
            
            % Indivdual run projections
            temp_vol = cell(expt.Nruns, 1);  
            tempRunProj = cell(expt.Nruns, 2, projParam.Nz);  
            runBinLims = cell(1,expt.Nruns);
            runProjCell = cell(expt.Nruns, 3, projParam.Nz); 
            runVolCell = cell(expt.Nruns, 3); 
            catCell = cell(projParam.Nz, 3);
            projParam.path.run.(projParam.sbx_type{t}).z = runProjCell; 
            projParam.path.run.(projParam.sbx_type{t}).vol = runVolCell;
            projParam.path.cat.(projParam.sbx_type{t}).z = runProjCell; 
            projParam.path.cat.(projParam.sbx_type{t}).vol = cell(1,3);

            for runs = 1:expt.Nruns
                tempName = sprintf('%s_run%i',expt.name, runs);
                if expt.Nplane == 1
                    [~, tempRunProj(runs,:,1), runBinLims{runs}, tempRawProjPath] = WriteSbxPlaneTif(expt.sbx.(projParam.sbx_type{t}), catInfo, 1, 'chan','both', 'firstScan',expt.scanLims(runs)+projParam.bin+1, 'Nscan',expt.Nscan(runs)-projParam.bin, ...
                        'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'sbxType',projParam.sbx_type{t}, 'RGB',false, 'dir',runsDir, 'name',tempName, 'overwrite',projParam.overwrite ); % , 'sbxType',projParam.sbx_type{t}
                    projParam.path.run.raw.z(runs,:,1) = tempRawProjPath;
                else
                    
                    % Z projections
                    if exist(expt.sbx.(projParam.sbx_type{t}), 'file')
                        [tempProj, runBinLims{runs}, tempProjPath] = WriteSbxZproj(expt.sbx.(projParam.sbx_type{t}), catInfo, 'z',projParam.z, 'chan','both', 'dir',runsDir, 'name',tempName, ...
                            'sbxType',projParam.sbx_type{t}, 'projType',projParam.type, 'monochrome',true, 'RGB',true, 'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), ...
                            'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite); % expt.Nscan(runs)
                        
                        % Reshape results to tempRunProj form
                        if iscell(tempProj)
                            for Z = 1:projParam.Nz
                                if size(tempProj{Z}, 4) == 1 && strcmpi(projParam.color{1}, 'green')
                                    tempRunProj(runs,2,:) = tempProj;
                                    projParam.path.run.(projParam.sbx_type{t}).z(runs,2,:) = tempProjPath;
                                elseif size(tempProj{Z}, 4) == 1 && strcmpi(projParam.color{1}, 'red')
                                    tempRunProj(runs,1,:) = tempProj;
                                    projParam.path.run.(projParam.sbx_type{t}).z(runs,2,:) = tempProjPath;
                                else
                                    tempRunProj{runs,1,Z} = tempProj{Z}(:,:,:,1);
                                    tempRunProj{runs,2,Z} = tempProj{Z}(:,:,:,2);
                                    projParam.path.run.(projParam.sbx_type{t}).z(runs,:,:) = permute(tempProjPath, [3,2,1]);
                                end
                            end
                        else % WriteSbxZproj returns tempProj as an array if only one set of z is provided
                            if size(tempProj, 4) == 1 && strcmpi(projParam.color{1}, 'green')
                                tempRunProj{runs,2,1} = tempProj;
                                projParam.path.run.(projParam.sbx_type{t}).z(runs,2,1) = tempProjPath;
                            elseif size(tempProj, 4) == 1 && strcmpi(projParam.color{1}, 'red')
                                tempRunProj{runs,1,1} = tempProj;
                                projParam.path.run.(projParam.sbx_type{t}).z(runs,1,1) = tempProjPath;
                            else
                                tempRunProj{runs,1,1} = tempProj(:,:,:,1);
                                tempRunProj{runs,2,1} = tempProj(:,:,:,2);
                                projParam.path.run.(projParam.sbx_type{t}).z(runs,:,1) = permute(tempProjPath, [3,2,1]);
                            end
                        end
                    end
                    % Volume tifs
                    if projParam.vol
                        if exist(expt.sbx.(projParam.sbx_type{t}), 'file')
                            [temp_vol{runs}, ~, projParam.path.run.(projParam.sbx_type{t}).vol(runs,:)] = WriteSbxVolumeTif(expt.sbx.(projParam.sbx_type{t}), catInfo, 'chan','both', 'dir',runsDir, 'name',tempName, 'type',projParam.sbx_type{t}, 'monochrome',true, 'RGB',true,...
                                'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'binT',projParam.bin, 'overwrite',projParam.overwrite);
                        end
                    end
                end

                % Account for temporal binning in time vectors
                if projParam.bin == 1
                    projParam.Tproj(runs) = Tscan(runs);
                else
                    for b = flip(1:size(runBinLims{runs},1))
                        projParam.Tproj{runs}(b) =  mean( Tcat(runBinLims{runs}(b,1):runBinLims{runs}(b,2)) );
                    end
                end
            end

            % Concatenate runs
            if expt.Nruns > 1
                fprintf('\nGenerating concatenated projections')
                if expt.Nplane == 1
                    for chan = find(any(~cellfun(@isempty, tempRunProj(:,:,1))))
                        if exist(expt.sbx.(projParam.sbx_type{t}), 'file')
                            projParam.path.cat.(projParam.sbx_type{t}).z{chan,1} = sprintf('%s%s_%s_%s.tif', projParam.dir, expt.name, projParam.sbx_type{t}, chanName{chan} );
                            if ~exist(projParam.path.cat.(projParam.sbx_type{t}).z{chan,1}, 'file')
                                WriteTiff(cat(3, tempRunProj{:,chan,1}), projParam.path.cat.(projParam.sbx_type{t}).z{chan,1} );
                            end
                        end
                    end
                else

                    % Concatenated z-projections
                    for Z = 1:projParam.Nz
                        for chan = find(all(~cellfun(@isempty, tempRunProj(:,:,Z))))
                            projParam.path.cat.(projParam.sbx_type{t}).z{Z,chan} = sprintf('%s%s_%s_z%i-%i_%s_%s.tif', projParam.dir, expt.name, projParam.type, projParam.z{Z}(1), projParam.z{Z}(end), projParam.sbx_type{t}, chanName{chan} );
                            if ~exist(projParam.path.cat.(projParam.sbx_type{t}).z{Z,chan}, 'file')
                                WriteTiff(cat(3, tempRunProj{:,chan,Z}), projParam.path.cat.(projParam.sbx_type{t}).z{Z,chan} );
                            end
                        end
                    end
                    
                    % Concatenated Volume tifs
                    if projParam.vol
                        if exist(expt.sbx.(projParam.sbx_type{t}), 'file')
                            projParam.path.cat.(projParam.sbx_type{t}).vol{3} = sprintf('%s%s_%s_vol_RGB.tif', projParam.dir, expt.name, projParam.sbx_type{t} );
                            if ~exist(projParam.path.cat.(projParam.sbx_type{t}).vol{3},'file') || projParam.overwrite
                                tempCat = uint16(cat(5, temp_vol{:}));
                                for chan = flip(1:size(tempCat,4))
                                    projParam.path.cat.(projParam.sbx_type{t}).vol{chan} = sprintf('%s%s_%s_vol_%s.tif', projParam.dir, expt.name, projParam.sbx_type{t}, chanName{chan} );
                                    if ~exist(projParam.path.cat.(projParam.sbx_type{t}).vol{chan},'file') || projParam.overwrite
                                        fprintf('\nWriting %s', projParam.path.cat.(projParam.sbx_type{t}).vol{chan})
                                        bfsave( tempCat(:,:,:,chan,:), projParam.path.cat.(projParam.sbx_type{t}).vol{chan});
                                    end

                                    chanLower = prctile(reshape(tempCat(:,:,:,chan,:),[],1), 5);
                                    chanUpper = max(reshape(tempCat(:,:,:,chan,:),[],1)); %prctile(stackChan{chan}(:), 1);
                                    %fprintf('\nRescaling %s channel : [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper);
                                    tempCat(:,:,:,chan,:) = rescale(tempCat(:,:,:,chan,:), 0, 2^8-1, 'inputMin',chanLower, 'inputMax',chanUpper);
                                end
                                fprintf('\nWriting %s', projParam.path.cat.(projParam.sbx_type{t}).vol{3})
                                bfsave( uint8(tempCat), projParam.path.cat.(projParam.sbx_type{t}).vol{3});
                            end
                        end
                    end
                end
            end
        else
            fprintf('\n%s does not exist - skipping', expt.sbx.(projParam.sbx_type{t}))
        end
    end
    projParam.Nbin = cellfun(@numel, projParam.Tproj);
    projParam.totBin = sum(projParam.Nbin);
    projParam.binLims = [0,cumsum(projParam.Nbin)];

    % Save the parameters, including Tproj
    fprintf('\nSaving %s', projParamPath)
    save(projParamPath, 'projParam');
else
    fprintf('\nLoading %s', projParamPath)
    load(projParamPath, 'projParam');
end