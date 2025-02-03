clc; clear all;

% TODO --- Set path to the SBX file to be z-projected
input_sbx = 'V:\2photon\Simone\Simone_Macrophages\Pf4Ai9-2\006\Pf4Ai9-2_240509_006.sbx';

% TODO --- set which planes to project
zProj = {3:6}; %{2:6, 19:23};

% TODO --- set max or mean projection (default = 'mean')
projType = 'mean'; %'mean'

% TODO --- set which channels you want ('both' will try to project both separately if 2 channels exist)
projChan = 'both'; %'both'

[projData, binLims, projTifPath] = WriteSbxZproj(input_sbx, [], 'sbxType','raw', 'z',zProj, 'chan',projChan, 'binT',1); %'sbxType','raw',  'Nscan',100

%{
, 'dir',runsDir, 'name',tempName, 'sbxType','raw', 'projType',projParam.type, 'monochrome',true,...
    'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite
%}