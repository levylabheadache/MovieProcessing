function runInfo = GetRunInfo(expt, dataDir)

for r = flip(1:expt.Nruns)
    runInfo(r) = MakeInfoStruct( dataDir, expt.mouse, expt.date, expt.runs(r), expt.fov );
end

end

