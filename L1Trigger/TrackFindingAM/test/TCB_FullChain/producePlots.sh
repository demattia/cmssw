cp /home/demattia/fdata_demattia/CMSSW_AMPR_Files/L1TrackNtuplePlot.C .
mkdir TrkPlots
root -l -b -q L1TrackNtuplePlot.C+\(\"Ntuple\", 13\)
cp /home/demattia/fdata_demattia/CMSSW_AMPR_Files/RunFake.C .
mkdir FakePlots
root -l -b -q RunFake.C+\(\"Ntuple\"\)
