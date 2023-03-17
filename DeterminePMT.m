function [usePMTind, usePMTname, useChan] = DeterminePMT(requestColor, sbxInfo)
% Which PMTs are available (PMT gain > 0)?
PMTname = {'green','red'};
pmtGain = [sbxInfo.config.pmt0_gain, sbxInfo.config.pmt1_gain];
pmtOn = pmtGain > 0; 
pmtAvail = PMTname(pmtOn);
% Which PMTs were requested?
if ~isempty(pmtAvail)
    if strcmpi(requestColor, 'both')
        pmtRequest = {'green','red'};
        if sbxInfo.nchan == 1, pmtRequest = pmtRequest{1}; end
    else
        pmtRequest = requestColor;
    end
    usePMTind = find(strcmpi(pmtRequest, pmtAvail));
    if isempty(usePMTind) 
        warning('Requested PMT (%s) not available, using %s instead!', requestColor, pmtAvail{1});
        usePMTind = find(strcmpi(pmtAvail{1}, PMTname));
    end
else
    error('No PMTs available!');
end
if numel(usePMTind) == 1
    usePMTname = PMTname{usePMTind};
else
    usePMTname = 'both';
end
chanInd = [2,1];
useChan = chanInd(usePMTind);

end

