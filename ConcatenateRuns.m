function ConcatenateRuns(expt, runInfo, catInfo, varargin )
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'runInfo', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
addParameter( IP, 'sbx', 'sbxz', @ischar )
addParameter( IP, 'chunkSize', 100, @isnumeric )
addParameter( IP, 'ext', 'sbxcat', @ischar )
addParameter( IP, 'overwrite', false, @islogical ) % for scanbox, 1 = green, 2 = red. -1 = both
parse( IP, expt, runInfo, catInfo, varargin{:} ); % mouse, exptDate,
sbxType = IP.Results.sbx;
catExt = IP.Results.ext; % 
chunkSize = IP.Results.chunkSize;
overwrite = IP.Results.overwrite;
catName = expt.name; % sprintf('%s_FOV%i', , expt.fov);
catPathRoot = sprintf('%s%s', expt.dir, catName);
catSbxPath = strcat(catPathRoot, '.', catExt);
if expt.Nruns > 1
    sbxPath = cell(1,expt.Nruns);
    if ~exist(catSbxPath,'file') || overwrite
        [~, ~, useChan] = DeterminePMT('both', catInfo); % which PMTs were used?
        Nchan = numel(useChan);
        % Write the sbxcat file
        fprintf('\n     Writing %s\n', catSbxPath); tic
        rw = SbxWriter(catSbxPath, catInfo, catExt, true); 
        w = waitbar(0, sprintf('writing %s',catExt));
        for runs = 1:expt.Nruns %expt.runs
            % Load the run stack, in chunks
            [runChunkLims, NrunChunk, runChunkLength] = MakeChunkLims(1, runInfo(runs).Nscan, runInfo(runs).Nscan, 'allowPartial',true, 'size',chunkSize);
            sbxPath{runs} = sprintf('%s%s.%s', runInfo(runs).dir, runInfo(runs).fileName, sbxType ); % sbx_affine  sbxz
            fprintf('\n   Loading %s (%i chunks at a time)... ', sbxPath{runs}, chunkSize); tic
            if expt.Nplane > 1
                for chunk = 1:NrunChunk
                    runChunkStack = readSBX(sbxPath{runs}, runInfo(runs), runChunkLims(chunk,1), runChunkLength(chunk), -1, []); % [c,x,y,z,t]
                    if Nchan == 2
                        runChunkStack = reshape(runChunkStack, [size(runChunkStack,[1,2,3]), prod(size(runChunkStack,[4,5]))]); % rw expects this form
                    elseif Nchan == 1
                        runChunkStack = reshape(runChunkStack, [size(runChunkStack,[1,2]), prod(size(runChunkStack,[3,4]))] );
                    end
                    rw.write( runChunkStack ); % Write the chunk to sbxcat
                end
            else
                for chunk = 1:NrunChunk
                    runChunkStack = readSBX(sbxPath{runs}, runInfo(runs), runChunkLims(chunk,1), runChunkLength(chunk), -1, []); % [c,x,y,t]
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
    % If the experiment consists of only one run, just copy it to the main folder
    [~, source_path] = FileFinder(runInfo.dir, 'contains',sbxType);
    if ~isempty(source_path)
        fprintf('\nSingle run experiment: copying %s to %s', source_path{1}, catSbxPath);
        copyfile(source_path{1}, catSbxPath);
    else
        error('\nNo %s file found in %s', sbxType, catSbxPath)
    end
end
end