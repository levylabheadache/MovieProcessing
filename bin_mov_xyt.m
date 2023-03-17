function binned = bin_mov_xyt(mov, factor, factor_t, removelines)
% binned = bin_mov(mov, factor) spatially bins the Y-by-X-by-nFrames movie
% by an integer factor: Each new pixel is the sum of factor-by-factor old
% pixels. If the dimensions of mov do not match, pixels are discarded at 
% the edges.
%modified by MA to do temporal binning as well

[y, x, Nframes, ~] = size(mov);

warning('on');
if mod(x, factor)
% 	warning('Pixels in X-dimension are not an integer multiple of the binning factor -- discarding %d Y-columns.', mod(x, factor));
	x = x - mod(x,factor);
end

if mod(y, factor)
% 	warning('Pixels in Y-dimension are not an integer multiple of the binning factor -- discarding %d Y-rows.', mod(y, factor));
	y = y - mod(y,factor);
end


if mod(Nframes, factor_t)
% 	warning('Pixels in T-dimension are not an integer multiple of the binning factor -- discarding %d T-frames.', mod(Nframes, factor_t));
	Nframes = Nframes - mod(Nframes,factor_t);
end

mov = mov(1:y,1:x,1:Nframes);

[y, x, z] = size(mov);

% Turn movie into 2-D vector
mov = reshape(mov, y, x*z);

[m,n] = size(mov);

% Bin along columns:
mov = sum(reshape(mov,factor,[]) ,1);

% Bin along rows:
mov = reshape(mov,m/factor,[]).'; %Note transpose
mov = sum(reshape(mov,factor,[]) ,1);
mov = reshape(mov,n/factor,[]).'; %Note transpose 

% Turn back into original shape:
binned = reshape(mov, y/factor, x/factor, z);

[y, x, z] = size(binned);
mov = reshape(binned,y*x,z)';
mov = sum(reshape(mov,factor_t,[]),1);
%binned = reshape(mov,y,x,z/factor_t);
binned = shiftdim(reshape(mov,z/factor_t,y,x),1);



if removelines == 1
    %now do highpass filtering
    tmp = squeeze(mean(binned,3));
    sz2 = size(tmp);
    tmpA = tmp(round(sz2(1)/8):round(sz2(1)*7/8),:);
    tmpB = tmp(:,round(sz2(1)/8):round(sz2(1)*7/8));
    
    
    %tmpA = tmp;
    %tmpB = tmp;
    
    tmp2 = (tmp - ones(sz2(1),1)*median(tmpA,1) + median(median(tmpA,1)) - median(tmpB,2)*ones(1,sz2(2)) + median(median(tmpB,2)));
    %figure; imagesc(tmp2,[mean(mean(tmp2)) + std(nonzeros(tmp2))*[-1 1]]);
    
    for count = 1:size(binned,3)
        tmp3 = squeeze(binned(:,:,count));
        binned(:,:,count) = (tmp3 - ones(sz2(1),1)*median(tmpA,1) + median(median(tmpA,1)) - median(tmpB,2)*ones(1,sz2(2)) + median(median(tmpB,2)));
    end
end

%tmp3 = (tmp - ones(sz(1),1)*mean(tmp,1) + mean(mean(tmp,1)) - mean(tmp,2)*ones(1,sz(2)) + mean(mean(tmp,2)));
%imagesc(tmp3,[mean(mean(tmp3)) + std(nonzeros(tmp3))*[-1 1]]);

 

clear mov;





