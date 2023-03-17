function [ref_run, ref_scan] = DetermineReference(expt, Tscan, loco, eligible_runs, minStillDur) % refFrame
% Find the longest epoch of stillness within the set of eligbile runs, and use that period as the reference point for registration
if isempty(eligible_runs), eligible_runs = expt.runs; end
%minStillDur = 30; 
trim = 5; 
%minTrim = 0.5;
%Nepoch = 0;

% Try progressively shorter minimum epoch durations until we find at least one epoch
ep_longest = [];
while isempty(ep_longest)%Nepoch == 0
    fprintf('\nLooking for a reference epoch: min dur = %2.1f s, trim %2.1f s from each side', minStillDur, trim)
    close all;
    [stillEpoch, ~, stillSumm] = BinStillEpochs(expt, Tscan, loco, [], [], [], 'criterion','speed', 'show',false, 'minStillDur',minStillDur, 'trim',trim); % true
    %Nepoch =  stillSumm.Nepoch;
    
    % Find the longest epoch within the eligible runs
    %{
    longest_dur = 0;
    for ep = find(ismember(stillSumm.epoch_run, eligible_runs)) % 1:stillSumm.Nepoch
        if stillSumm.dur_full(ep) > longest_dur
            ep_longest = ep;
            longest_dur = stillSumm.dur_full(ep);
        end
    end
    %}
    % which run had the longest still epoch?
    tempDur = stillSumm.dur_full;
    tempDur(~ismember(stillSumm.epoch_run, eligible_runs)) = NaN; % suppress ineligible epochs
    [~,ep_longest] = max(tempDur);

    if ~isempty(ep_longest)
        ref_run = stillSumm.epoch_run(ep_longest);
        [~,refRunLongestEp] = max(stillEpoch(ref_run).dur_full);
        ref_scan = stillEpoch(ref_run).scan_trim{refRunLongestEp};
        %{
        ref_frame = stillEpoch(ref_run).scan_trim{refRunLongestEp};
        % Convert frames to scans, if necessary
        if expt.Nplane > 1
            ref_scan = unique(ceil((ref_frame)/expt.Nplane));
            ref_scan([1,end]) = []; % remove edge
        else
            ref_scan = ref_frame;
        end
        %}
    end
    % Lower the threshold for the next iteration, if needed  
    minStillDur = minStillDur-1;
    if trim >= minStillDur/4, trim = minStillDur/10; end 
    %trim = max([trim-0.5, minTrim]);
end
%{
if stillSumm.Nepoch > 0
    if ~isnan(expt.csd)
        [~,longestEpRunInd] = max(stillSumm.dur_trim(stillSumm.epoch_run < expt.csd)); % find the longest pre-CSD epoch
    else
        [~,longestEpRunInd] = max(stillSumm.dur_trim); % find the longest pre-CSD epoch
    end
    refRun = stillSumm.epoch_run(longestEpRunInd);
    [~,longestEpInd] = max(stillEpoch(refRun).Nscan_trim);
    refScan = stillEpoch(refRun).scan_trim{longestEpInd};
end
%}
end