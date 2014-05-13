#!/bin/tcsh -f
#
# Script to batch merge all files of a directory. 
# Expects a run directory as input.
#

# Avoid wild cards
set noglob

##### CONFIGURATION ##############################################
# Where to merge the files to (will be copied over at the end)
set TOPWORKDIR=/scratch/`whoami`

set dir=$1
set run=$2
set numfiles="-0"
if ( "$3" != "" ) then
  set numfiles="$3"
endif
set scratch=$TOPWORKDIR/$run
set merged=$TOPWORKDIR/$run.root
set remote=`echo $dir | sed 's@\/*$@@'`.root # Remove trailing slashes...
set srmPath="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/"
##################################################################


############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

##### MONITORING/DEBUG INFORMATION ###############################
set DATE_START=`date +%s`
echo "Job started at " `date`


##### RUNNING ####################################################

### Check arguments
if ( "$1" == "" ) then
   echo "Usage: $0 <run directory>"
   exit -1
endif


# Set the CMSSW environment (for ROOT, mainly...)
source $VO_CMS_SW_DIR/cmsset_default.csh
setenv SCRAM_ARCH slc5_amd64_gcc462 
cmsenv
# Fix problem with DCache and ROOT
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH

echo Processing $numfiles files from $dir

# Get list of files in directory
set files = `lcg-ls --vo cms --srm-timeout=6000 --connect-timeout=6000 -T srmv2 --nobdii $srmPath$dir | awk '{print "dcap://t3se01.psi.ch:22125"$1}' | head --lines=$numfiles | tr '\n' ' '`
echo Merging to $merged from:
echo $files

# Copy files locally
set lfiles = ()
mkdir -pv $scratch
foreach file ( $files )
   set lfile = $scratch/$file:t
   dccp "$file" $lfile
   set lfiles = ($lfiles $lfile)
end

# Now merge
hadd -f $merged $lfiles
if ( $? != 0 ) then
  echo Merging failed
  rm -vf $merged $lfiles
  exit -1
endif

# Move merged file to DCache (remove if already exists)
echo "Removing local files"
rm -vf $lfiles
ls $scratch
rmdir $scratch
echo "Removing remote file (if existing) $remote"
srmrm -debug=true $srmPath$remote
echo "Copying merged file to SE"
lcg-cp -b -D srmv2 -U srmv2 $merged $srmPath$remote && rm -vf $merged

echo ALL DONE



###########################################################################
set DATE_END=`date +%s`
@ RUNTIME = $DATE_END - $DATE_START
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0
