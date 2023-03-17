clc; clear all;

% TODO --- Set path to the SBX file to be z-projected
input_sbx = 'D:\2photon\Anna\NaVAi9AG\002\NaVAi9AG_221122_002.sbxfix';

% TODO --- set which planes to project
zProj = {2:6}; %{2:6, 19:23};

% TODO --- set max or mean projection (default = 'mean')
projType = 'mean'; %'mean'

% TODO --- set which channels you want ('both' will try to project both separately if 2 channels exist)
projChan = 'green'; %'both'

[projData, binLims, projTifPath] = WriteSbxZproj(input_sbx, [], 'sbxType','raw', 'z',zProj, 'chan',projChan, 'binT',1); %'sbxType','raw',  'Nscan',100

%{
, 'dir',runsDir, 'name',tempName, 'sbxType','raw', 'projType',projParam.type, 'monochrome',true,...
    'firstScan',expt.scanLims(runs)+1, 'Nscan', expt.Nscan(runs), 'edge',projParam.edge, 'scale',projParam.scaleFactor, 'binT',projParam.bin, 'overwrite',projParam.overwrite
%}