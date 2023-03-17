function CorrectData3D(sbxPath, sbxInfo, shiftPath, varargin) % , tforms_optotune , Nchunks
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'sbxInfo', @isstruct )
addRequired(IP, 'shiftPath', @ischar )
addOptional(IP, 'refChan', 'green', @ischar ) % for scanbox, PMT1 = green, 2 = red. -1 = both
addParameter(IP, 'chunkSize', 15, @isnumeric)
%addParameter(IP, 'type', 'median', @ischar); %options are 'median' or 'mean' for projecting the reference volume
addParameter(IP, 'blur', 1, @isnumeric); %width of gaussian to blur for DFT reg
addParameter(IP, 'keep', 0.95, @isnumeric); %what proportion of frame to accountf or with shifts
addParameter(IP, 'anchor',0, @isnumeric); %
addParameter(IP, 'scale', 2, @isnumeric ) % 
addParameter(IP, 'edges',[0,0,0,0], @isnumeric); % [left, right, top, bottom]
addParameter(IP, 'save', true, @islogical);
parse(IP, sbxPath, sbxInfo, shiftPath, varargin{:}); % tforms_optotune, 
refChan = IP.Results.refChan;
chunkSize = IP.Results.chunkSize;
scaleFactor = IP.Results.scale;
blurFactor = IP.Results.blur;
keepFactor = IP.Results.keep;
saveToggle = IP.Results.save;
edges = IP.Results.edges;
%refType = IP.Results.type;
anchor = IP.Results.anchor;

[refPMT, ~] = DeterminePMT(refChan, sbxInfo); % PMT1 = green, PMT2 = red  refPMTname
Nx = sbxInfo.sz(1); Ny = sbxInfo.sz(2);
NxCrop = Nx - edges(3) - edges(4);
NyCrop = Ny - edges(1) - edges(2);
[chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.Nscan, 'size',chunkSize );

RS = cell(1,Nchunk);  CS = cell(1,Nchunk);  ZS = cell(1,Nchunk);
tic;
w = waitbar(0,'DFT registration'); %H = parfor_progressbar(Nchunk,'DFT registration'); % 
for c = 1:Nchunk
    % load chunk of optotune-corrected data
    raw_chunk = readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkLength(c), refPMT ); % readSBX(path, info, k, N, pmt, optolevel)
    raw_chunk = raw_chunk(edges(3)+1:end-edges(4), edges(1)+1:end-edges(2), :);
    raw_chunk = reshape(raw_chunk, NxCrop, NyCrop, sbxInfo.Nplane, []);
    raw_chunk = imresize(raw_chunk, 1/scaleFactor); 

    % rectify each volume with dft
    RS0 = zeros(sbxInfo.Nplane,chunkSize); CS0 = zeros(sbxInfo.Nplane,chunkSize); 
    if sbxInfo.Nplane > 1
        chunk_reg0 = nan(size(raw_chunk));
        parfor i = 1:chunkSize
            %[RS0(:,i), CS0(:,i), chunk_reg0(:,:,:,i)] = RectifyStackDFT( raw_chunk(:,:,:,i), round(sbxInfo.Nplane/2), 4); % DFT_rect( raw_chunk(:,:,:,i), , 4);
            [RS0(:,i), CS0(:,i), chunk_reg0(:,:,:,i)] = RectifyStack( raw_chunk(:,:,:,i), round(sbxInfo.Nplane/2)); % DFT_rect( raw_chunk(:,:,:,i), , 4);
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
    RS2 = shifts(:,1)*scaleFactor;
    CS2 = shifts(:,2)*scaleFactor;
    ZS1 = shifts(:,3);

    % combine Row and Column Shifts from DFT registrations above
    RS(:,c) = {RS0*scaleFactor + RS1*scaleFactor + repmat(RS2',[sbxInfo.Nplane,1])};
    CS(:,c) = {CS0*scaleFactor + CS1*scaleFactor + repmat(CS2',[sbxInfo.Nplane,1])};
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
RS_chunk = interchunk_shifts(:,1)*scaleFactor;
CS_chunk = interchunk_shifts(:,2)*scaleFactor;
ZS_chunk = interchunk_shifts(:,3);
RS = [RS{:}];  CS = [CS{:}];  ZS = [ZS{:}]; %convert local shift correction cells to matrix form

%stretch the intra-chunk corrections to apply to every frame
RS_chunk = imresize(RS_chunk',size(RS),'nearest');
CS_chunk = imresize(CS_chunk',size(CS),'nearest');
ZS_chunk = imresize(ZS_chunk',size(ZS),'nearest');

%save the DFT and optotune transformations
if saveToggle
    save(shiftPath, 'RS', 'CS', 'ZS', 'RS_chunk', 'CS_chunk', 'ZS_chunk', 'intermediate_shifts', 'scaleFactor', 'ref_all', ...
        'chunkSize', 'anchor', 'edges', 'scaleFactor', 'keepFactor', 'blurFactor', '-mat', '-v7.3'); % ,'tforms_optotune_full'
    %load( shiftPath )
end
    
end