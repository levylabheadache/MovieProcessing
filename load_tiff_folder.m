function mov = load_tiff_folder(fDir)
%load_tiff_folder(fDir) Loads all tiffs in fDir and concatenates them into
%one long movie.

% Get total number of tiffs and frames in folder:
nTiffs = 0;
nFrames = 0;
for f = dir(fDir)'
	[~, ~, ext] = fileparts(f.name);
	if strcmp(ext, '.tif')
		fInfo = imfinfo(fullfile(fDir, f.name)  ,'tiff');
		nFrames = nFrames + length(fInfo);
		nTiffs = nTiffs +1;
	end
end

nCols = fInfo(1).Width;
nRows = fInfo(1).Height;

% Pre-allocate RAM for speed:
% (This method is much faster than using ZEROS)
mov = [];
mov(nRows, nCols, nFrames) = 0;

% Load individual tiff files:
multiWaitbar('Loading TIFFs from folder', 0);
frame_start = 1;

for f = dir(fDir)'
	[~, ~, ext] = fileparts(f.name);
	if ~strcmp(ext, '.tif')
		continue
	end
	
	mov_part = load_tiff_nobar(fullfile(fDir, f.name));
	frame_end = (frame_start - 1) + size(mov_part,3);
	mov(:,:,frame_start:frame_end) = mov_part;
	
	frame_start = frame_end + 1;
	
	multiWaitbar('Loading TIFFs from folder', frame_end/nFrames);
end

multiWaitbar('Loading TIFFs from folder', 'close');