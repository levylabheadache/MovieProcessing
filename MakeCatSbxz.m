function MakeCatSbxz(sbxPath, sbxInfo)
%OBSOLETE?
% Check if the sbxz file already exists
[~, sbxzPath] = FileFinder( sbxInfo.dir, 'type','mat', 'type','sbxz' );
if isempty(sbxzPath)
    % Load interpolation results
    [~, zinterpPath] = FileFinder( sbxInfo.dir, 'type','mat', 'contains','zinterp' );
    if ~isempty(zinterpPath)
        fprintf('\nLoading %s', zinterpPath{1})
        load(zinterpPath{1},'-mat');
        % Generate .sbxz file
        rw = SbxWriter(sbxPath, sbxInfo, '.sbxz');
        w = waitbar(0,'Generating sbxz');
        for scan = 1:sbxInfo.totScan
            % Load one scan at a time
            data_scan = readSBX(sbxPath, sbxInfo, scan, 1, -1, []);
            % Apply Zshifts
            if sbxInfo.nchan == 1
                data_scan = ApplyZShiftInterpolateFBS(data_scan, ZS(:,scan), CS(:,scan), RS(:,scan)); %toc % apply shifts, this step is slow
                data_scan = reshape(data_scan, [size(data_scan,[1,2]), prod(size(data_scan,[3,4]))] );
            else
                data_scan = permute(data_scan, [2,3,4,5,1]); % ApplyZShiftInterpolateFBS expects spatial dimensions first
                for c = 1:2 % (parallel seems to be slower than serial)
                    data_scan(:,:,:,:,c) = ApplyZShiftInterpolateFBS(data_scan(:,:,:,:,c), ZS(:,scan), CS(:,scan), RS(:,scan));
                end
                data_scan = permute(data_scan, [5,1,2,3,4]); %toc % regwriter expects color dim first
                data_scan = reshape(data_scan, [size(data_scan,[1,2,3]), prod(size(data_scan,[4,5]))] );
            end
            % Write result to .sbxz
            rw.write(data_scan);
            waitbar(scan/sbxInfo.totScan, w);
        end
        toc
        delete(w);
        rw.delete;
    else
        error('\nNo zinterp file found!'); 
    end
else
    fprintf('\n%s already exists - returning\n', sbxzPath{1})
end
end