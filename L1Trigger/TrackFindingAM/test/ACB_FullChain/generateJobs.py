import os

from listFiles import list_files

list_files()

index = 0
for input_file_name in open("filelist.txt"):
    cfg_name = "AM_CB_"+str(index)+".py"
    output_cfg = open(cfg_name, "w")
    for line in open("AM_CB.py"):
        output_cfg.write(line.replace("INDEX", str(index)).replace("INPUT_FILE_NAME", input_file_name.rstrip('\n')))
    output_cfg.close()

    job_file_name = "jobFile_"+str(index)+".slrm"
    output_job_file = open(job_file_name, "w")
    output_job_file.write("""#!/bin/bash
#SBATCH -J PCA
#SBATCH -p hepx
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
cd """+os.getcwd()+"""
cmsRun """+cfg_name+"""

echo "Ended at `date` on `hostname`"

exit 0
""")
    output_job_file.close()

    index += 1

    print "sbatch "+job_file_name
    os.system("sbatch "+job_file_name)
