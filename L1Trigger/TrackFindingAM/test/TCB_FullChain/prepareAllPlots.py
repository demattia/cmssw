import os

def prepare_plots(file_type):
    print "Producing plots in "+file_type+"_ntuples"
    # print "cd "+file_type+"_ntuples"
    # print "hadd Ntuple.root Ntuple_*.root"
    # print "rm Ntuple_*.root"
    # print "rm L1TrackNtupleMaker*.py"
    # print "rm jobFileNtuple*"
    # print "cp /home/demattia/fdata_demattia/CMSSW_AMPR_Files/producePlots.sh ."
    # print "./producePlots.sh"
    # print "cd -"
    os.system("""cd """+file_type+"""_ntuples
pwd
hadd Ntuple.root Ntuple_*.root
rm Ntuple_*.root
rm L1TrackNtupleMaker*.py
rm jobFileNtuple*
cp /home/demattia/fdata_demattia/CMSSW_AMPR_Files/producePlots.sh .
./producePlots.sh
cd -""")


prepare_plots("TC")
# prepare_plots("CB")
prepare_plots("FIT")
prepare_plots("DR")
