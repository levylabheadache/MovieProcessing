function MakeSbxPartial(expt, sbxInfo, varargin) 
% Apply only the some components of the affine transform
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'sbxInfo', @isstruct )
addParameter( IP, 'outPath', '', @ischar )
addParameter( IP, 'sbx', true, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'suppress', {'scale', 'shear'}, @iscell)
parse( IP, expt, sbxInfo, varargin{:} ); 
sbxOutputPath = IP.Results.outPath;
suppress = IP.Results.suppress; 

% Load the affine transforms
tformPath = sprintf('%s%s_affine_tforms.mat', expt.dir, expt.name); % %strcat(fDir,'\',fName,'_affine_tforms.mat'); % '_regTforms.mat'
fprintf('\nLoading %s... ', tformPath );
load( tformPath ); % 'sbxPath', 'sbxInfo', 'refVol', 'params', 'tforms_
% Determine sbx file names
sbxInputPath = sprintf('%s%s.sbxz', expt.dir, expt.name);
if isempty(sbxOutputPath)
    sbxOutputPath = sprintf('%s%s.sbx_partial', expt.dir, expt.name); % sprintf('%s\\%s.sbx_partial', fDir, fName);
end
% Determine which transforms to suppress
suppressTrans = false; suppressScale = false; suppressShear = false;
if any(strcmpi(suppress, 'trans'))
    suppressTrans = true;
end
if any(strcmpi(suppress, 'scale'))
    suppressScale = true;
end
if any(strcmpi(suppress, 'shear'))
    suppressShear = true;
end
if (suppressTrans & suppressScale & suppressShear) || (~suppressTrans & ~suppressScale & ~suppressShear), error('No types suppressed!'); end

% Apply the transforms and generate sbx_partial
Nrow = sbxInfo.sz(1); 
Ncol = sbxInfo.sz(2); 
Nchan = sbxInfo.nchan;
Nscan = sbxInfo.totScan; 
Nplane = sbxInfo.Nplane;
tic
w = waitbar(0,'writing .sbx\_partial');
rw = pipe.io.RegWriter(sbxOutputPath, sbxInfo, '.sbx_partial', true);
if Nchan == 1
    fprintf('\n     Writing monochrome partially-registered sbx file'); 
    aff_chunk = zeros(Nrow, Ncol, Nplane);
    for s = 1:Nscan
        % Apply affine transformation to each plane of the scan
        input_chunk = readSBX(sbxInputPath, sbxInfo, Nplane*(s-1)+1, Nplane, 1, []); % pipe.imread(sbxInputPath, Nplane*(s-1)+1, Nplane, params.refChan, []);
        for z = 1:Nplane % parfor is actually slower
            tform_partial = tforms_all{z,s};
            if suppressTrans, tform_partial.T(3,[1,2]) = [0,0]; end % suppress translation
            if suppressScale, tform_partial.T([1,5]) = [1,1]; end % suppress scaling
            if suppressShear, tform_partial.T([2,4]) = [0,0]; end % suppress shearing
            aff_chunk(:,:,z) = imwarp(input_chunk(:,:,z), tform_partial, 'OutputView',imref2d([Nrow,Ncol]));
        end
        % Write the results to sbx_affine
        rw.write(squeeze(uint16(aff_chunk))); %rw.write(squeeze(uint16(tempScan)));
        waitbar(s/Nscan, w);
    end
else
    error('Only made for monochrome data currently');
end
rw.delete;
delete(w);
toc
end