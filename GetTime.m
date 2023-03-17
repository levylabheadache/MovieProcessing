function [Tscan, infoStruct, Tcat] = GetTime(infoStruct)
% Get trans-movie time vectors
Nruns = numel(infoStruct);
Tscan = cell(1,Nruns); dTtrans = zeros(1,Nruns);
for r = 1:Nruns
    % infoStruct(r).Nscan = floor([infoStruct(r).nframes]/infoStruct(r).otlevels);
    dT = infoStruct(r).Nplane/infoStruct(r).framerate; 
    dTtrans(r) = 3600*24*(datenum(infoStruct(r).timestamp) - datenum(infoStruct(1).timestamp));
    Tscan{r} = dT*(0:infoStruct(r).Nscan-1)' + dTtrans(r);
end
Tcat = vertcat(Tscan{:});

end

