function MakeSbxZ(sbxPath, sbxInfo, shiftPath, varargin)
% Load interpolation results
load(shiftPath,'-mat');
RS = RS1 + RS2 + RS3 + RS_chunk;
CS = CS1 + CS2 + CS3 + CS_chunk;
ZS = ZS1 + ZS_chunk;

tic
fprintf('\nLoading %s', sbxPath); tic;
sbx_data = readSBX(sbxPath, sbxInfo, 1, sbxInfo.Nscan, -1, []); toc

fprintf('\nApplying z interpolation shifts');  tic;
if sbxInfo.nchan == 1
    sbx_data = ApplyZShiftInterpolateFBS(sbx_data, ZS, CS, RS); toc % apply shifts, this step is slow
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2]), prod(size(sbx_data,[3,4]))] ); 
else
    sbx_data = permute(sbx_data, [2,3,4,5,1]); % ApplyZShiftInterpolateFBS expects spatial dimensions first
    % Interpolate each color, serial (parallel seems to be slower than serial)
    for c = 1:2
        sbx_data(:,:,:,:,c) = ApplyZShiftInterpolateFBS(sbx_data(:,:,:,:,c), ZS, CS, RS);
    end
    % Interpolate each color 
    sbx_data = permute(sbx_data, [5,1,2,3,4]); toc % regwriter expects color dim first
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2,3]), prod(size(sbx_data,[4,5]))] ); 
end
toc

fprintf('\nWriting .sbxz'); tic
rw = SbxWriter(sbxPath, sbxInfo, '.sbxz'); % rw = pipe.io.RegWriter(sbxPath, sbxInfo, '.sbxz', true); 
rw.write(sbx_data);
rw.delete;
toc
end