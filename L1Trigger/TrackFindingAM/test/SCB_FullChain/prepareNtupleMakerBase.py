import os
from os import listdir
from os.path import isfile

def prepare_ntuple_maker_cfg(cfg_type):
    files_path = os.getcwd()
    only_files = [f for f in listdir(files_path) if f.find("_output") != -1 and f.find(".root") != -1]
    
    template_cfg = open("L1TrackNtupleMaker"+cfg_type+"_template.py", "r")
    new_cfg = template_cfg.readlines()
    cfg_dir = cfg_type+"_ntuples"
    os.system("mkdir -p Backup")
    os.system("rm -rf Backup/"+cfg_dir)
    os.system("mv "+cfg_dir+" Backup/")
    os.system("mkdir -p "+cfg_dir)
    
    i=0
    for f in only_files:
        input_file = ("\"file:"+files_path+"/" + f + "\"\n")
    
        output_cfg = open(cfg_dir+"/"+"L1TrackNtupleMaker"+cfg_type+"_"+str(i)+".py", "w")
    
        for line in new_cfg:
            line = line.replace("INPUT_FILE_LIST", input_file)
            line = line.replace("INDEX", str(i))
            output_cfg.write(line)
    
        output_cfg.close()
    
        output_job_file_name = cfg_dir+"/"+"jobFileNtuple_"+str(i)
        output_job_file = open(output_job_file_name, "w")
        output_job_file.write("""#!/bin/bash
#SBATCH -J PCA
#SBATCH -p background-4g
#SBATCH -o job.txt
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8000mb
# # SBATCH --mem-per-cpu=2000mb
#SBATCH -n1

echo "Starting at `date` on `hostname`"

echo "SLURM_JOBID=$SLURM_JOBID"
echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAH_TASK_ID=$SLURM_ARRAY_TASK_ID"

cd /home/demattia/Seb/CMSSW_6_2_0_SLHC28_patch1/src
source /home/hepxadmin/cmssw/cmsset_default.sh
eval `scramv1 runtime -sh`
cd """+files_path+"/"+cfg_dir+"/"+"""
cmsRun L1TrackNtupleMaker"""+cfg_type+"_"+str(i)+""".py

echo "Ended at `date` on `hostname`"

exit 0
""")
        output_job_file.close()
        
        print "sbatch "+output_job_file_name
        os.system("sbatch "+output_job_file_name)
    
        i += 1
