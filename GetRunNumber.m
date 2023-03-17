function runNumber = GetRunNumber(x)
runNumber = str2double(x(strfind(x, 'run')+3:end));
end