function segParams = GetSegParams(sbxInfo)

paramsPath = sprintf('%s%s_seg_params.mat', sbxInfo.dir, sbxInfo.exptName ); 
if exist(paramsPath, 'file')
    fprintf('\nLoading %s', paramsPath);
    load(paramsPath, 'IP');
    segParams = IP.Results;
else
    error('%s does not exist', paramsPath);
end


end

