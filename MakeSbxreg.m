function MakeSbxreg(sbxInfo) % , sbxInputPath, sbxOutputPath, params

% Define input and output sbx paths
sbxInputPath = sbxInfo.path;
sbxOutputPath = sprintf('%s%s.sbxreg', sbxInfo.dir, sbxInfo.exptName);


% Load regTforms_corrected and params
[~, tformPath] = FileFinder(sbxInfo.dir, 'type','mat', 'contains','regTforms');
fprintf('\nLoading %s', tformPath{1}); tic;
tformStruct = load(tformPath{1});
if isfield(tformStruct, 'regTform_correct')
    fprintf('\nUsing corrected regTforms')
    regTform = tformStruct.regTform_correct;
else
    fprintf('\nUsing uncorrected regTforms')
    regTform = tformStruct.regTform;
end
params = tformStruct.params;
toc

% Apply the transforms and generate sbx_affine
[chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.totScan, 'N',10);
imRef = imref2d([sbxInfo.sz(1), sbxInfo.sz(2)]);
w = waitbar(0,'writing .sbxreg');
rw = SbxWriter(sbxOutputPath, sbxInfo, '.sbxreg', true);
fprintf('\n     Writing registered sbx file');
tic
if sbxInfo.nchan == 1
    [pmt, ~] = DeterminePMT(params.refChan, sbxInfo);
    tic
    for chunk = 1:Nchunk
        data_chunk = readSBX(sbxInputPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), pmt, []);
        if sbxInfo.Nplane == 1
            % Single plane, single color
            for s = 1:chunkLength(chunk)
                data_chunk(:,:,s) = imwarp(data_chunk(:,:,s), regTform{1,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
            end
        else
            % Multi plane, single color
            for s = 1:chunkLength(chunk)
                for z = 1:sbxInfo.Nplane
                    data_chunk(:,:,z,s) = imwarp(data_chunk(:,:,z,s), regTform{z,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                end
            end
            data_chunk = reshape(data_chunk, [size(data_chunk,[1,2]), prod(size(data_chunk,[3,4]))]);
        end

        rw.write(data_chunk);
        waitbar(chunk/Nchunk, w);
        toc
    end
else
    tic
    for chunk = 1:Nchunk
        data_chunk = readSBX(sbxInputPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk), -1, []);
        if sbxInfo.Nplane == 1
            % Single plane, multi color
            data_chunk = permute(data_chunk, [2,3,4,1]);
            for s = 1:chunkLength(chunk)
                for chan = 1:2
                    data_chunk(:,:,s,chan) = imwarp(data_chunk(:,:,s,chan), regTform{1,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                end
            end
            data_chunk = permute(data_chunk, [4,1,2,3]);
        else
            % multi plane, multi color
            data_chunk = permute(data_chunk, [2,3,4,5,1]);
            for s = 1:chunkLength(chunk)
                for chan = 1:2
                    for z = 1:sbxInfo.Nplane
                        data_chunk(:,:,z,s,chan) = imwarp(data_chunk(:,:,z,s,chan), regTform{z,chunkLims(chunk,1)+s-1}, 'OutputView',imRef);
                    end
                end
            end
            data_chunk = permute(data_chunk, [5,1,2,3,4]);
            data_chunk = reshape(data_chunk, [size(data_chunk,[1,2,3]), prod(size(data_chunk,[4,5]))]);
        end
        rw.write(data_chunk);
        waitbar(chunk/Nchunk, w);
        toc
    end
end
rw.delete;
delete(w);
toc
end