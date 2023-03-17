%% Remove bad imaging data - in case the 2P recording had any issue half way through that affected unpacked movie

targetDir = 'D:\2photon\AB39\220922_FOV1\014\';
targetPath = 'D:\2photon\AB39\220922_FOV1\014\AB39_220922_014.sbx';
tempInfo = MakeInfoStruct( targetPath );
WriteSbxPlaneTif(targetPath, tempInfo, 1, 'dir',targetDir, 'name','AB39_220922_014', 'verbose',true, 'chan','green', 'binT',1, 'overwrite',true);
tempInfo = FixSBX(targetPath, tempInfo, 'flip',flipZ, 'proj',true, 'overwrite',false, 'scans',1:1312);