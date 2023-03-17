function x = readSBX(path, info, firstScan, Nscan, pmt, z)
%
if nargin < 2
    [pathDir, pathName, ~] = fileparts(path); % pathExt
    [~,infoPath] = FileFinder(pathDir, 'type','mat', 'criteria',@(x)(strcmp(x,pathName))); % 
    fprintf('\nNo info structure was provided. Loading %s', infoPath{1})
    info = MakeInfoStruct( infoPath{1} );
end

% Which scan to begin reading from?
if nargin < 3, firstScan = 1; end  

% How many scans, including the first, to read?
if nargin < 4 || Nscan <= 0,  Nscan = info.totScan;  end
if firstScan + Nscan - 1 > info.totScan, Nscan = info.totScan - firstScan + 1; end % Make sure that we don't search beyond the end of the file
if rem(Nscan,1) ~= 0, error('Nscan must be a positive integer'); end

% Determine which PMT(s) to read
if nargin < 5,  pmt = -1; end % 1 = green, 2 = red, -1 = all
if pmt == -1 
    if info.nchan == 2
        pmt = [1,2]; 
    elseif info.nchan == 1  
        pmt = find([info.config.pmt0_gain, info.config.pmt1_gain], 1, 'first'); % if PMT set to -1, but only one channel was recorded, default to green
    end 
end

% If optotune was used, which planes to read?
if nargin < 6,  z = []; end
if ~isempty(z) && z ~= 1 && ~info.optotune_used
    error('Optotune was not used for this file.');
end

% Open the file
fileID = fopen( path );
if fileID == -1, error(['Cannot read file ' path]); end

if isempty(z)
    % Read all planes from the desired scans
    seekStatus = fseek(fileID, (firstScan-1)*info.Nplane*info.nsamples, 'bof'); % fseek(info.fid, k*info.nsamples, 'bof');
    if seekStatus == 0
        x = fread(fileID, Nscan*info.Nplane*info.nsamples/2, 'uint16=>uint16');
        x = reshape(x, [info.nchan, info.sz(2), info.recordsPerBuffer, info.Nplane, Nscan]);
    else
        error(ferror(fileID));
    end
else
    % Read only the selected plane z
    bufWidth = info.nchan*info.sz(2)*info.recordsPerBuffer;
    x = zeros(bufWidth, Nscan, 'uint16'); 
    k = 0;
    for s = firstScan:firstScan+Nscan-1
        k = k+1;
        seekStatus = fseek(fileID, (info.Nplane*(s-1)+(z-1))*info.nsamples, 'bof'); 
        if seekStatus == 0
            x(:,k) = fread(fileID, info.nsamples/2, 'uint16=>uint16'); %info.fid % x(n*bufWidth+1:(n+1)*bufWidth)
            %{
            im = intmax('uint16') - permute( reshape(x(:,k), [info.nchan, info.sz(2), info.recordsPerBuffer]), [3,2,1] );
            ch1 = im(:,:,1);  subplot(1,2,1); imshow(ch1, []);
            ch2 = im(:,:,2);  subplot(1,2,2); imshow(ch2, []);
            %}
        end
    end
    %x = x(:); 
    x = reshape(x, [info.nchan, info.sz(2), info.recordsPerBuffer, Nscan]);
end
fclose(fileID);
x = intmax('uint16') - permute(x, [1 3 2 4 5]);
x = squeeze(x(pmt,:,:,:,:)); % squeeze the output if a single PMT is called

end
