function correctedVolume = ApplyZShiftInterpolateFBS(Volume, ZShifts, ColumnShifts, RowShifts)

%   APPLYXYZSHITINTERPOLATE: apply XY shifts per plane with a parfor loop
%   using a cubic interpolation
%
%   Inputs:
%     Volume -- 4D matrix of uint16 or other, dim (x,y,z,t)
%     ZShifts -- 2D matrix of doubles, dim (z,t)
%     ColumnShifts -- 2D matrix of doubles, dim (z,t)
%     RowShifts -- 2D matrix of doubles, dim (z,t)
%   Outputs:
%     correctedVolume -- 4D matrix of uint16 or other, dim (x,y,z,t)

correctedVolume = zeros(size((Volume)));
WidthCorr = floor((size(Volume, 3) - size(ZShifts, 1)));

for t = 1:size(Volume,4)
    fri = 1:size(Volume,3) ;   % plane serie to be inteprolated
    MiniVol = double(Volume(:,:,:,t));
    for plane = 1:size(Volume,3)
        if plane > WidthCorr && plane <= size(Volume,3)-WidthCorr
            C = ColumnShifts(plane,t);
            R = RowShifts(plane,t);
            slice = correctedVolume(:,:,plane,t);
            correctedVolume(:,:,plane,t)=imtranslate(slice,[C, R]);
            fri(plane) = plane + ZShifts(plane,t);
        else
            correctedVolume(:,:,plane,t) = double(Volume(:,:,plane,t));
        end
    end
    % interpolation
    ZSerieNonReg = timeseries(MiniVol, fri); % zseries non reg volume
    ZSerieReg = resample(ZSerieNonReg, 1:size(Volume, 3));
    % resample volume with linear interpolation on every integer plane
    correctedVolume(:,:,:,t) = ZSerieReg.Data;
end
end
