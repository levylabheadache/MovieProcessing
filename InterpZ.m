function InterpZ(sbxPath, sbxInfo, varargin)
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar )
addRequired(IP, 'sbxInfo', @isstruct )
addOptional(IP, 'shiftpath', '', @ischar )
addParameter(IP, 'chunkSize', 15, @isnumeric )
addParameter(IP, 'refChan', 'green', @ischar ) % for scanbox, 1 = green, 2 = red. -1 = both
addParameter(IP, 'refType', 'mean', @ischar); %options are 'median' or 'mean' for projecting the reference volume
addParameter(IP, 'scale', 4, @isnumeric ) %
addParameter(IP, 'edges',[0,0,0,0], @isnumeric); % [left, right, top, bottom]
addParameter(IP, 'planescorr', 3, @isnumeric); %how many planes up and down to search for Z DFT reg
addParameter(IP, 'blur', 1, @isnumeric); %width of gaussian to blur for DFT reg
addParameter(IP, 'keep',0.95, @isnumeric); %what proportion of frame to accountf or with shifts
%addParameter(IP, 'save',true, @islogical); %save the DFT and optotune shifts
%addParameter(IP, 'debug',false, @islogical); %save intermediate tiffs
parse( IP, sbxPath, sbxInfo, varargin{:} );
shiftPath = IP.Results.shiftpath;
chunkSize = IP.Results.chunkSize;
refChan = IP.Results.refChan;
[refPMT, ~] = DeterminePMT(refChan, sbxInfo); % PMT1 = green, PMT2 = red
refType = IP.Results.refType;
edges = IP.Results.edges;
scale = IP.Results.scale;
planescorr = IP.Results.planescorr;
blurFactor = IP.Results.blur;
keepingFactor = IP.Results.keep;
%saveToggle = IP.Results.save;
%debugToggle = IP.Results.debug;

fDir = fileparts(sbxPath);
Nx = sbxInfo.sz(1); % - edges(3) - edges(4);
Ny = sbxInfo.sz(2); % - edges(1) - edges(2);
[chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.Nscan, 'size',chunkSize );
RS = cell(3,Nchunk); CS = cell(3,Nchunk); ZS = cell(1,Nchunk);
tic
w = waitbar(0,'Z interpolation...');
for c = 1:Nchunk %parfor
    % Load current chunk
    raw_chunk = readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkLength(c), refPMT ); % readSBX(sbxPath, sbxInfo, chunkFrames*(chunk-1)+1, chunkFrames, refChan, [] );  
    raw_chunk = reshape(raw_chunk, Nx, Ny, sbxInfo.Nplane, []);
    raw_chunk = raw_chunk(edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:,:);
    raw_chunk = imresize(raw_chunk,1/scale);

    %if debugToggle, pipe.io.writeTiff(raw_chunk, strcat(fDir,filesep,'raw_chunk.tif') ); end

    % First round of XY DFT registration
    %make single target volume for entire chunk  %ref1 = defineReference(raw_chunk, size(raw_chunk,4), refType);
    if strcmpi(refType, 'mean')
        ref1 = mean(raw_chunk, 4);
    elseif strcmpi(refType, 'median')
        ref1 = median(raw_chunk, 4);
    end
    [RS1, CS1] = DetermineXYShiftsFBS(raw_chunk, blurFactor, keepingFactor, ref1); %calculate DFT shift between each volume and target vol
    reg_chunk1 = ApplyXYShiftsFBS(raw_chunk,RS1,CS1); %apply the shift (need it for subsequent steps)
    %if debugToggle,  pipe.io.writeTiff(reg_chunk1, strcat(fDir,filesep,'reg_chunk1.tif'));  end
    raw_chunk = []; %clear raw chunk to save space

    % First round of Z interpolation registration
    % make new reference volume based on XY DFT registration % ref2 = defineReference(reg_chunk1, size(reg_chunk1,4), refType); 
    if strcmpi(refType, 'mean')
        ref2 = mean(reg_chunk1, 4);
    elseif strcmpi(refType, 'median')
        ref2 = median(reg_chunk1, 4);
    end
    [RS2, CS2, ZS1] = InterpolateZshift(ref2, reg_chunk1, planescorr);
    %[RS2,CS2,ZS1] = ComputeZshiftInterpolateFBS(ref2, reg_chunk1, planescorr, edges); %calculate interpolated zshifts by fitting polynomial
    reg_chunk2 = ApplyZShiftInterpolateFBS(reg_chunk1, ZS1, CS2, RS2); %apply zshift
    %if debugToggle, pipe.io.writeTiff(reg_chunk2,strcat(fDir,filesep,'reg_chunk2.tif'));  end
    %find NaNs and replace with zeros
    nan_idx = find(isnan(reg_chunk2));
    reg_chunk2(nan_idx) = 0;
    %clear registration from step 1
    reg_chunk1 = [];

    % Final round of XY DFT registration
    %make new reference volume based on XYZ DFT (4)  %ref3 = defineReference(reg_chunk2, size(reg_chunk2,4), refType);
    if strcmpi(refType, 'mean')
        ref3 = mean(reg_chunk2, 4);
    elseif strcmpi(refType, 'median')
        ref3 = median(reg_chunk2, 4);
    end
    [RS3, CS3] = DetermineXYShiftsFBS(reg_chunk2, blurFactor, keepingFactor, ref3); %calculate DFT shift between each volume and target vol
    %{
    if debugToggle
        reg_chunk3 = ApplyXYShiftsFBS(reg_chunk2,RS3,CS3);
        pipe.io.writeTiff(reg_chunk3,strcat(fDir,filesep,'reg_chunk3.tif'));
    end
    %}
    reg_chunk2 = []; %clear to save space

    % Combine Row and Column Shifts from previous steps
    RS(:,c) = {RS1,RS2,RS3};
    CS(:,c) = {CS1,CS2,CS3};
    ZS(c) = {ZS1};

    % Retain final reference volumes for stitching later
    ref_all{c} = ref3;

    waitbar(c/Nchunk, w) 
end
delete(w);
toc
%Fix intra-chunk discontinuities
ref_final = ref_all{1};
%concatenate all of the reference volumes from above
ref_cat = ref_all{1};
for i = 2:Nchunk
    ref_cat = cat(4,ref_cat,ref_all{i});
end

%Compute XYZ shifts from each chunk's reference volume
[RS_chunk,CS_chunk,ZS_chunk] = InterpolateZshift(ref_final, ref_cat, planescorr); %ComputeZshiftInterpolateFBS(ref_final, ref_cat, planescorr, edges);
clearvars ref_cat %ref_cat = [];

%convert local shift correction cells to matrix form
CS1 = [CS{1,:}];
CS2 = [CS{2,:}];
CS3 = [CS{3,:}];
RS1 = [RS{1,:}];
RS2 = [RS{2,:}];
RS3 = [RS{3,:}];
ZS1 = [ZS{:}];

%stretch the intra-chunk corrections to apply to every frame
RS_chunk = imresize(RS_chunk,size(RS1),'nearest');
CS_chunk = imresize(CS_chunk,size(CS1),'nearest');
ZS_chunk = imresize(ZS_chunk,size(ZS1),'nearest');
%{
if debugToggle
    tic
    raw_chunk = pipe.imread(sbxPath,chunkFrames*(1-1)+1, chunkFrames,refChan,[]);
    ch = ApplyXYShiftsFBS(raw_chunk,RS1(:)*scale,CS1(:)*scale);
    fprintf('one shift applied\n');
    ch = ApplyZShiftInterpolateFBS(ch,ZS1(:),CS2(:)*scale,RS2(:)*scale);
    fprintf('two shifts applied\n');
    ch = ApplyXYShiftsFBS(ch,RS3(:)*scale,CS3(:)*scale);
    fprintf('three shifts applied\n');
    ch = ApplyXYShiftsFBS(ch,RS_chunk(:)*scale,CS_chunk(:)*scale);
    fprintf('four shifts applied\n');
    toc
    pipe.io.writeTiff(ch,strcat(fDir,filesep,'chunk_reg_piecewise.tif'));
    fprintf('written final DFT\n');
end
%}

%save the DFT and optotune transformations
if ~isempty(shiftPath) %saveToggle
    save(shiftPath,'RS1','CS1','ZS1','RS2','CS2','RS3','CS3','RS_chunk','CS_chunk','ZS_chunk','scale', '-mat');
end

end