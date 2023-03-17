function downsampled = downsampleWithAvg(mov, fDown)
% downsampleWithAvg reduces the number of frames in movie MOV by the factor
% FDOWN, by averaging blocks of successive frames.

if fDown < 0 || fDown - round(fDown) ~= 0 || fDown > size(mov,3)
	error('Downsampling factor must be a positive integer between 1 and the number of frames in mov.');
end

% Truncate movie such that the number of frames is evenly divisible by
% fDown
mov = mov(:,:,1:end-mod(size(mov,3),fDown));

% Split movie array into separate blocks according to fDown to facilitate
% averaging
mov_d = zeros(size(mov,1),...
				size(mov,2),...
				size(mov,3)/fDown,...
				fDown);
		
for i = 1:fDown
	mov_d(:,:,:,i) = mov(:,:,1:fDown:end);
end

% Average along fourth dimension
downsampled = mean(mov_d, 4);