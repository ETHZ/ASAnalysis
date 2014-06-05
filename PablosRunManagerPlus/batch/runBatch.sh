#!/bin/tcsh -f
#
# PhysQC processing script
#
# Expects output directory and full input file path as input. 
#

# Avoid wild cards
set noglob

##### CONFIGURATION ##############################################
# Where to store local results (will be copied over at the end)
set TOPWORKDIR=/scratch/`whoami`
##################################################################


############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N p_run

# Run time soft and hard limits hh:mm:ss
# soft=CPU time; hard=Wallclock time
#$ -l s_rt=40:00:00,h_rt=42:00:00

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

##### MONITORING/DEBUG INFORMATION ###############################
set DATE_START=`date +%s`
echo "Job started at " `date`


##### RUNNING ####################################################

### Check arguments
if ( $# < 2 ) then
   echo "Usage: $0 <exe> <output directory> <file1> [... <filen>]"
   exit -1
endif

set exe = $1
shift
set dir = $1
shift
set scratch = $TOPWORKDIR/$dir
set files = "$*"

# Set the CMSSW environment (for ROOT, mainly...)
source $VO_CMS_SW_DIR/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc434
cmsenv
# Fix problem with DCache and ROOT
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH

echo Processing $dir from $files with $exe

./$exe -o $scratch $files

# Remove output directory, if it exists
rm -rvf ./$dir

# Move output directory to local directory, and merged file to DCache
mv -v $scratch .


###########################################################################
set DATE_END=`date +%s`
@ RUNTIME = $DATE_END - $DATE_START
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0
