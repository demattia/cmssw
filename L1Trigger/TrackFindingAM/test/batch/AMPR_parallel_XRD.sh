#!/bin/bash

###########################################
#
# Main script for parallel CMSSW AM-based pattern 
# recognition in batch on a multi-core machine
# 
# !!!Working on EDM files!!!
#
# !!!Requires the installation of GNU parallel executable on your machine!!!
#    !!!!!!! http://www.gnu.org/software/parallel/ !!!!!!!
#
# If you cannot use parallel, set p6 to 1 in the launch command
# this could also apply when you have one file and one bank only
#
# The jobs themselves are launched by PR_processor_parallel.sh
#
# source AMPR.sh p1 p2 p3 p4 p5 p6 p7 p8
# with:
# p1 : The directory containing the data file you want to analyze (best is to copy them beforehand on the machine scratch area)
# p2 : Name of the SE subdirectory where you will store the data
# p3 : The directory where you will retrieve the bank files, the pattern reco will
#      run over all the pbk files contained in this directory
# p4 : How many events per input data file? 
# p5 : How many events per job (should be below p3...)?
# p6 : The global tag name
# p7 : How many cores you want to use in parallel (if one then parallel is not used)
# p8 : How many events per job to process
#
# For more details, and examples, have a look at:
# 
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP V.2)
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 01/03/2016
#
# Script tested with release CMSSW_6_2_0_SLHC27
#
###########################################


# Here we retrieve the main parameters for the job 

INPUTDIR=${1}   # Directory where the input root files are
OUTPUTDIR=${2}  # Name of the output directory in the SE
BANKDIR=${3}    # Directory where the bank (.pbk) files are
NTOT=${4}       # How many events per data file?
NPFILE=${5}     # How many events per job?
GTAG=${6}       # Global tag
NCORES=${7}     # #cores
NFILES=${8}     # #files per job

###################################
#
# The list of parameters you can modify is here
#
###################################

export LFC_HOST=lyogrid06.in2p3.fr

# I/O Directory info on the grid
BAREPATH=/dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC
SRMPATH=srm://$LFC_HOST/$BAREPATH
XRDPATH=root://$LFC_HOST/$BAREPATH

# The parallel command

parallel=/gridsoft/ipnls/parallel/bin/parallel

###########################################################
###########################################################
# You are not supposed to touch the rest of the script !!!!
###########################################################
###########################################################


# The SE directory containing the input EDM files

INDIR_GRID=$SRMPATH/$INPUTDIR
INDIR_XROOT=$XRDPATH/$INPUTDIR

# The SE directory containing the output EDM file with the PR output
OUTDIR_GRID=$SRMPATH/$OUTPUTDIR
OUTDIR_XROOT=$XRDPATH/$OUTPUTDIR

MATTER=`basename $OUTPUTDIR`

OUTDIRTMP=/data/viret/AMProc/TMPDAT_$MATTER
OUTDIRSTO=/data/viret/AMProc/TMP_$MATTER

echo 'The data will be read from directory: '$INDIR_GRID
echo 'The final pattern reco output files will be written in: '$OUTDIR_GRID

# Following lines suppose that you have a certificate installed on lxplus. To do that follow the instructions given here:
#
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideLcgAccess
#

voms-proxy-init --voms cms --valid 100:00 -out $HOME/.globus/gridproxy.cert
export X509_USER_PROXY=${HOME}/.globus/gridproxy.cert

lfc-mkdir -p $BAREPATH/$OUTPUTDIR

# Then we setup some environment variables

PACKDIR=$PWD            # This is where the scripts are created 
RELEASEDIR=$CMSSW_BASE  # This is where the release is installed

mkdir $OUTDIRTMP
mkdir $OUTDIRSTO

# We loop over the data directory in order to find all the files to process

ninput=0	 
nsj=0
npsc=$NFILES

echo "#\!/bin/bash" > global_stuff_${MATTER}.sh

for ll in `lcg-ls $INDIR_GRID | grep root` 
do   
    l=`basename $ll`
    echo $l
    i=0
    j=$NPFILE

    val=`expr $ninput % $npsc`

    if [ $val = 0 ]; then

	nsj=$(( $nsj + 1))

	echo "source $PACKDIR/run_${nsj}_${MATTER}.sh"  >> global_stuff_${MATTER}.sh

	echo "#\!/bin/bash" > run_${nsj}_${MATTER}.sh

	if [ $NCORES = 1 ]; then
	    echo "$PACKDIR/run_PRMERGE_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "$PACKDIR/run_FMERGE_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "$PACKDIR/run_FIT_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "$PACKDIR/run_RM_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	else
	    echo "${parallel} -j ${NCORES} < $PACKDIR/run_PRMERGE_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "${parallel} -j ${NCORES} < $PACKDIR/run_FMERGE_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "${parallel} -j ${NCORES} < $PACKDIR/run_FIT_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "$PACKDIR/run_RM_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	fi

	echo "#\!/bin/bash" > run_PRMERGE_${nsj}_${MATTER}.sh
	echo "#\!/bin/bash" > run_FMERGE_${nsj}_${MATTER}.sh
	echo "#\!/bin/bash" > run_FIT_${nsj}_${MATTER}.sh
	echo "#\!/bin/bash" > run_RM_${nsj}_${MATTER}.sh

	chmod 755 run_${nsj}_${MATTER}.sh
	chmod 755 run_PRMERGE_${nsj}_${MATTER}.sh
	chmod 755 run_FMERGE_${nsj}_${MATTER}.sh
	chmod 755 run_FIT_${nsj}_${MATTER}.sh
	chmod 755 run_RM_${nsj}_${MATTER}.sh

    fi

    ninput=$(( $ninput + 1))

    echo 'Working with file '$l

    # First look if the file has been processed

    OUTM=`echo $l | cut -d. -f1`

    OUTF=${OUTM}"_with_AMPR.root"
    OUTE=${OUTM}"_with_FIT.root"
    OUTD=${OUTM}"_extr.root"

    processed=0
    section=0

    AMPR_FILE=${OUTDIR_GRID}/$OUTF
    FIT_FILE=${OUTDIR_GRID}/$OUTE
    EXTR_FILE=${OUTDIR_GRID}/$OUTD

    dealF=`lcg-ls $AMPR_FILE | wc -l`
    dealE=`lcg-ls $FIT_FILE | wc -l`
    dealD=`lcg-ls $EXTR_FILE | wc -l`

    if [ $dealF != "0" ] && [ $dealE != "0" ] && [ $dealD != "0" ]; then
	
	echo "File "$l" has already processed, skip it..."
	continue;
    fi

    while [ $i -lt $NTOT ]
    do

	sec=0
        secdone=0
        section=$(( $section + 1))

	echo "$PACKDIR/PR_processor_parallel_p.sh PRMERGE ${i}_${j}.root $INDIR_XROOT $l ${OUTM}_ $RELEASEDIR $GTAG ${OUTM}_${i}_${j} $OUTDIRSTO $BANKDIR $i $NPFILE $OUTDIRTMP" >> run_PRMERGE_${nsj}_${MATTER}.sh

	i=$(( $i + $NPFILE ))
	j=$(( $j + $NPFILE ))

    done # End of loop over one input file
	
    #
    # Third step, all the merged files for the given input
    # file have been processed. Then launch the final merging 
    # 

    echo "$PACKDIR/PR_processor_parallel_p.sh  FINAL MERGED_${OUTM}_ $OUTDIRTMP $OUTDIRTMP $OUTF $RELEASEDIR ${OUTM} $OUTDIRSTO"  >> run_FMERGE_${nsj}_${MATTER}.sh
    echo "$PACKDIR/PR_processor_parallel.sh  FIT $OUTDIRTMP/${OUTF} $OUTE $OUTD $NTOT $OUTDIR_GRID $RELEASEDIR $GTAG ${OUTM} $OUTDIRSTO" >> run_FIT_${nsj}_${MATTER}.sh
    echo "cd $OUTDIRTMP" >> run_RM_${nsj}_${MATTER}.sh
    echo "rm *${OUTM}_*root" >> run_RM_${nsj}_${MATTER}.sh
done

chmod 755 global_stuff_${MATTER}.sh

