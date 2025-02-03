%% UnpackSBX extracts all sbx files within a folder - it is the FIRST step after you finish recording on the 2photon

% Clear any previous variables in the Workspace and Command Window to start fresh
clear; clc; close all; 

% TODO --- Set the directory of your 'sbx' files
mainDir = 'D:\2photon\Simone\Simone_Macrophages\Pf4Ai162-23\'; 
% add '\' to mainDir if it is not included
if mainDir(end) ~= '\'
    mainDir(end+1) = '\';
end
[sbxName, sbxPath] = FileFinder(mainDir, 'type','sbx');

% Flip Z - recording is performed backwards if optotune(3D) is used (bottom to top), so we should flip it
flipZ = true;

% Create S for unpacking backwards (from last run to first - setting up total number of runs facilitates and speeds up processing)
for s = flip(1:numel(sbxPath)) %1:numel(sbxPath) % flip(1:numel(sbxPath))
    sbxInfoPath = sbxPath{s};
    sbxInfoPath(end-2:end) = 'mat';
    
    % Move all files to a subfolder
    tempInfo = MakeInfoStruct( sbxInfoPath );

    % create subfolder name (00N)
    subDir = strcat(mainDir, tempInfo.fileName(end-2:end),'\');
    mkdir( subDir );
    [~,tempFiles] = FileFinder(mainDir, 'contains', tempInfo.fileName  ); %  , 'type', '*' , 'type','sbx'
    
    % move/copy files to subfolder
    for f = 1:numel(tempFiles)
        try
            movefile(tempFiles{f}, subDir, 'f') % try to move the file
        catch
            fprintf('\nmovefile failed, copying instead')
            copyfile(tempFiles{f}, subDir) 
        end
    end
    
    % Unpack the sbx file in the new subfolder - this function is for initial/basic analysis of the recordings
    [~,  newSBXpath] = FileFinder(subDir, 'type','sbx', 'keepExt',true);
    if ~tempInfo.optotune_used  %for 2D recording
        WriteSbxPlaneTif(newSBXpath{1}, tempInfo, 1, 'dir',subDir, 'name',sbxName{s}, 'verbose',true, 'chan','both', 'binT',8, 'overwrite',false);  % tempStack =
        WriteSbxProjection(newSBXpath{1}, tempInfo, 'verbose',true, 'chan','both', 'monochrome',true, 'RGB',true, 'overwrite',false); % projPath, 
    else  %for 3D recording
        tempInfo = FixSBX(newSBXpath{1}, tempInfo, 'flip',flipZ, 'proj',true, 'overwrite',false);
        fixPath = strcat(newSBXpath{1},'fix');   
        zProj = round(prctile(1:tempInfo.Nplane, 25)):round(prctile(1:tempInfo.Nplane, 75));
        WriteSbxZproj(fixPath, tempInfo, 'z',zProj, 'chan','both', 'overwrite',false, 'sbxType','raw', 'projType','mean', 'RGB',false, 'monochrome',true); % zprojPath,        
    end
end

% TODO --- Function movefile is not working on the server.
% 1 - Make sure you have a copy of these 3 files (.sbx .mat. .quad) saved on the \R:\Levylab\2photon
% 2 - Delete server versions that are duplicated
