%% ProjectionDemo creates tiff files with other frame range(s) - useful when the original frame selection does not include the structures that you need
%NOTE: THIS SCRIPT IS USEFUL FOR UNPACKING SINGLE FILES QUICKLY. USE GenerateExptProjections to write out projections for multiple runs from an experiment

% NOTE -- this script works only for 3D recordings
% If you need to binT a 2D recording, do ????

clear; clc; close all; % clear any previous variables in the Workspace and Command Window to start fresh

% Set path of .sbxfix previously generated with UnpackSBX
sbxPath = 'D:\2photon\Simone\Pf4Ai162-1\221122_FOV1\Pf4Ai162-1_221122_FOV1_run1\Pf4Ai162-1_221122_001.sbx';
%'D:\2photon\Simone\Pf4Ai162-1\221116_FOV3\Pf4Ai162-1_221116_FOV3_run1\Pf4Ai162-1_221116_003.sbxfix'; % exact path of the file to be projected
infoPath = []; % the sbxInfo structure, the path to the sbxInfo .mat file, or leave blank to determine automatically 
% Make the projection with default parameters (note, default doesn't write RGB projection. default params are estbalished in WriteSbxZproj function)
try
    WriteSbxZproj(sbxPath, infoPath); % volumetric data case
catch
    WriteSbxPlaneTif(sbxPath, infoPath, 1); % single plane case
end

% Make the projection with other parameters
zProj = {2:5}; % which set(s) of planes to project over
crop_edges = [60,60,20,20]; % crop this many pixels from the left/right,top/bottom edges
proj_name = 'Pf4Ai162-1_221116_003_crop_Tbin5'; % leave blank to set automatically; 'Pf4Ai162-1_221106_FOV2_run1_crop_Tbin5_XYbin4_scans2-102'
binT = 5; % temporal averaging - average every 5 scans together
%binXY = 4; % spatial averaging - average every square 4x4 pixels together
first_scan = 2; % first scan of the projection
Nscan = -1; % # of scans to project, set to -1 to use all, OR remove the field and value altogether from the function call; 100

% Project over selected sets of planes from volumetric data
WriteSbxZproj(sbxPath, infoPath, 'z',zProj, 'firstScan',first_scan, 'edges',crop_edges, 'chan','both', 'projType','mean', ...
    'name',proj_name, 'RGB',true, 'monochrome',true, 'overwrite',false); % zprojPath, tempInfo, 'binT',binT, 'scale',binXY,'Nscan',Nscan, 

% Project single-plane data
WriteSbxPlaneTif(sbxPath, infoPath, 1, 'firstScan',first_scan, 'edges',crop_edges, 'chan','both', 'projType','mean', ...
    'name',proj_name, 'RGB',true, 'monochrome',true, 'overwrite',false);