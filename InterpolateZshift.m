function[RowShifts,ColumnShifts,ZShifts] = InterpolateZshift(ReferenceVolumes, Volume, WidthCorr)
%   COMPUTEZSHIFTINTERPOLATE: calculate Z shifts
%
%   Inputs:
%     ReferenceVolumes -- 4D matrix of uint16 or other, dim (x,y,z/n,t)
%     Volume -- 2D matrix of doubles, dim (z,t)
%     WidthCorr -- 2D matrix of doubles, dim (z,t)
%     Edges -- list of int [xmin, xmax, ymin, ymax], Careful with the
%       resersed order using optotune
%   Outputs:
%     ColumnShifts -- 2D matrix of doubles, dim (z,t)
%     RowShifts -- 2D matrix of doubles, dim (z,t)
%     ZShifts -- 2D matrix of doubles, dim (z,t)

% Parameters
Size = size(Volume);
Size(4) = size(Volume,4);
nbChunck = size(ReferenceVolumes, 4); % nb chuncks of the reference
Chunck = floor(Size(4)/nbChunck); % nb frames by chunck of ref

% Downsample for speed and memory
Volume = Volume(1:2:end,1:2:end,:,:);
ReferenceVolumes = ReferenceVolumes(1:2:end,1:2:end,:,:);

% Preallocating space for output variables
RowShifts = zeros(Size(3),Size(4));
ColumnShifts = zeros(Size(3),Size(4));
ZShifts = zeros(Size(3),Size(4));

for t=1:Size(4) % for each time frame
    % pick corresponding reference
    reft = ReferenceVolumes(:,:,:,ceil(t/Chunck));
    % preallocate output vectors
    row_shift = zeros((Size(3)-WidthCorr),1);
    column_shift = zeros((Size(3)-WidthCorr),1);
    z_shift = zeros((Size(3)-WidthCorr),1);

    for j = WidthCorr+1:Size(3)-WidthCorr % j = considered plane
        Corrj = ones(1,2*WidthCorr+1)*NaN;

        for i = j-WidthCorr:j+WidthCorr % i = reference plane
            ref_slice = reft(:,:,i);
            slice = Volume(:,:,j,t);
            Corrj(1,i-j+WidthCorr+1) = corr2(ref_slice,slice);
        end
        [~,J] = max(Corrj);
        % Set interpolation vectors and degree regarding matrix size
        if  (J-5 > 0) & (J+5 <= size(Corrj))
            idx = 5;
        else
            idx = min([length(Corrj)-J, J-1]);
        end
        % Find Z shift based on spatial correlations
        x = J-idx:0.01:J+idx;
        FitOrder = idx;
        P = polyfit(J-idx:J+idx, Corrj(J-idx:J+idx),FitOrder);
        CorrelationFit = polyval(P, x);
        [~,I] = max(CorrelationFit); % max of the polynomial curve
        % recalcutate x and y shifts
        output = dftregistrationAlex(fft2(reft(:,:,j)),...
            fft2(Volume(:,:,j,t)),100);
        row_shift(j,1) = output(1);
        column_shift(j,1) = output(2);
        z_shift(j,1) = x(I)-WidthCorr-1;
    end
    % add shifts into output matrix
    RowShifts(:,t) = cat(1,row_shift, zeros(WidthCorr,1));
    ColumnShifts(:,t) = cat(1,column_shift, zeros(WidthCorr,1));

    % ensuring strict monotony, necessary?
    psteps = ones(length(z_shift),1).*(1:length(z_shift))';
    zaux = -z_shift + psteps;
    count = 0;
    while ~issorted(zaux) & count < 10
        count = count + 1;
        for plane = 2:size(z_shift)
            if zaux(plane) - zaux(plane-1) <= 0
                if ismember(count, [1,2,5,6,9,10]) % mod(count,2)==1;
                    zaux(plane-1) = NaN;
                else
                    zaux(plane) = NaN;
                end
            end
        end
        zaux = naninterp(zaux);
    end
    z_shift = zaux - psteps;
    ZShifts(:,t) = cat(1,-z_shift, zeros(WidthCorr,1));
end
end
