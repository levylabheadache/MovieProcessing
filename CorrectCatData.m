function CorrectExpt(expt, sbxInfo, regParam, varargin) % , tforms_optotune , Nchunks shiftPath,
IP = inputParser;
%addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'expt', @isstruct )
addRequired(IP, 'sbxInfo', @isstruct )
%addRequired(IP, 'shiftPath', @ischar )
%addOptional(IP, 'refChan', 'green', @ischar ) % for scanbox, PMT1 = green, 2 = red. -1 = both
addParameter(IP, 'chunkSize', 15, @isnumeric)
%addParameter(IP, 'type', 'median', @ischar); %options are 'median' or 'mean' for projecting the reference volume
addParameter(IP, 'blur', 1, @isnumeric); %width of gaussian to blur for DFT reg
addParameter(IP, 'keep', 0.95, @isnumeric); %what proportion of frame to accountf or with shifts
addParameter(IP, 'anchor',0, @isnumeric); %
addParameter(IP, 'save', true, @islogical);
parse(IP, sbxPath, sbxInfo, regParam, varargin{:}); % tforms_optotune,  shiftPath, 
chunkSize = IP.Results.chunkSize;
anchor = IP.Results.anchor;
blurFactor = IP.Results.blur;
keepFactor = IP.Results.keep;

[refPMT, ~] = DeterminePMT(regParam.refChan, sbxInfo); % PMT1 = green, PMT2 = red  refPMTname
Nx = sbxInfo.sz(1); Ny = sbxInfo.sz(2);
NxCrop = Nx-regParam.edges(3)-regParam.edges(4);
NyCrop = Ny-regParam.edges(1)-regParam.edges(2);
resize_factor = 1/regParam.binXY;
rectify_start = round(sbxInfo.Nplane/2);
[chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.Nscan, 'size',chunkSize );
RS = cell(1,Nchunk);  CS = cell(1,Nchunk);  ZS = cell(1,Nchunk);
tic;
w = waitbar(0,'DFT registration'); %H = parfor_progressbar(Nchunk,'DFT registration'); % 
for c = 1:Nchunk
    % load chunk of optotune-corrected data
    raw_chunk = readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkLength(c), refPMT ); % readSBX(path, info, k, N, pmt, optolevel)
    raw_chunk = raw_chunk(regParam.edges(3)+1:end-regParam.edges(4), regParam.edges(1)+1:end-regParam.edges(2), :);
    raw_chunk = reshape(raw_chunk, NxCrop, NyCrop, sbxInfo.Nplane, []);
    raw_chunk = imresize(raw_chunk, resize_factor); 

    % rectify each volume with dft
    RS0 = zeros(sbxInfo.Nplane,chunkSize); CS0 = zeros(sbxInfo.Nplane,chunkSize); 
    if sbxInfo.Nplane > 1
        chunk_reg0 = nan(size(raw_chunk));
        parfor i = 1:chunkSize
            [RS0(:,i), CS0(:,i), chunk_reg0(:,:,:,i)] = RectifyStack( raw_chunk(:,:,:,i), rectify_start) ; % DFT_rect( raw_chunk(:,:,:,i), , 4);
        end
    else
        chunk_reg0 = raw_chunk;
    end

    % first round of XY DFT registration
    ref1 = mean(chunk_reg0, 4); %defineReference(chunk_reg0, chunkSize, refType); % make single reference volume for entire chunk
    [RS1, CS1] = DetermineXYShiftsFBS(chunk_reg0, blurFactor, keepFactor, ref1); %calculate DFT shift between each volume and target vol
    chunk_reg1 = ApplyXYShiftsFBS(chunk_reg0,RS1,CS1); %apply the shift (need it for subsequent steps)

    % 3D DFT registration
    ref2 = mean(chunk_reg1, 4); %defineReference(chunk_reg1, chunkSize, refType);
    shifts = zeros(chunkSize,3);
    if sbxInfo.Nplane > 1
        parfor j = 1:chunkSize
            shifts(j,:) = dftregistration3D(fftn(ref2), fftn(chunk_reg1(:,:,:,j)), 2); % dftregistration3D(fftn(ref2), fftn(vol), 2);
        end
    end
    RS2 = shifts(:,1)*regParam.binXY;
    CS2 = shifts(:,2)*regParam.binXY;
    ZS1 = shifts(:,3);

    % combine Row and Column Shifts from DFT registrations above
    RS(:,c) = {RS0*regParam.binXY + RS1*regParam.binXY + repmat(RS2',[sbxInfo.Nplane,1])};
    CS(:,c) = {CS0*regParam.binXY + CS1*regParam.binXY + repmat(CS2',[sbxInfo.Nplane,1])};
    ZS(:,c) = {ZS1'};

    % save reference files for stitching later
    ref_all{c} = ref2;
    RS0_all{c} = RS0;
    CS0_all{c} = CS0;
    RS1_all{c} = RS1;
    CS1_all{c} = CS1;
    RS2_all{c} = RS2;
    CS2_all{c} = CS2;

    waitbar(c/Nchunk, w); % H.iterate(1); % 
end
close(w); % close(H); % 
toc
 
intermediate_shifts.RS0_all = [RS0_all{:}];
intermediate_shifts.RS1_all = [RS1_all{:}];
intermediate_shifts.RS2_all = [RS2_all{:}];
intermediate_shifts.CS0_all = [CS0_all{:}];
intermediate_shifts.CS1_all = [CS1_all{:}];
intermediate_shifts.CS2_all = [CS2_all{:}]; 

%Fix intra-chunk discontinuities
if anchor == 0, anchor = round(Nchunk/2); end
ref_final = ref_all{anchor};
interchunk_shifts = zeros(Nchunk,3);
if sbxInfo.Nplane > 1
    for j = 1:Nchunk
        interchunk_shifts(j,:) = dftregistration3D(fftn(ref_final),fftn(ref_all{j}),2);
    end
end
RS_chunk = interchunk_shifts(:,1)*regParam.binXY;
CS_chunk = interchunk_shifts(:,2)*regParam.binXY;
ZS_chunk = interchunk_shifts(:,3);
RS = [RS{:}];  CS = [CS{:}];  ZS = [ZS{:}]; %convert local shift correction cells to matrix form

%stretch the intra-chunk corrections to apply to every frame
RS_chunk = imresize(RS_chunk',size(RS),'nearest');
CS_chunk = imresize(CS_chunk',size(CS),'nearest');
ZS_chunk = imresize(ZS_chunk',size(ZS),'nearest');

% save the result transformations, and other relevant info
shiftPath = strcat(sbxInfo.dir, sbxInfo.exptName, '_dftshifts.mat');
fprintf('\nSaving %s', shiftPath)
save(shiftPath, 'regParam', 'chunkSize', 'anchor', 'keepFactor', 'blurFactor', 'RS', 'CS', 'ZS', 'RS_chunk', 'CS_chunk', 'ZS_chunk', 'intermediate_shifts', 'ref_all', '-mat', '-v7.3'); 

    
end