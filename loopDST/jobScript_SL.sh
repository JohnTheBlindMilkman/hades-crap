#!/bin/bash

# This is the generic jobscript to run jobs on GridEngine
#
# The script supports up to 7 parameters.
#
# The user specific part of the script is indicated below.

par1="no"
par2="no"
par3="no"
par4="no"
par5="no"
par6="no"
par7="no"

if [ $# -lt 1 ]
then
        echo "Usage : jobScript.sh par1 [par2] [par2] [par3] [par4] [par5] [par6] [par7]"
        sleep 3
        exit
fi

case "$#" in

1)  par1=$1
    ;;
2)  par1=$1
    par2=$2
    ;;
3)  par1=$1
    par2=$2
    par3=$3
    ;;
4)  par1=$1
    par2=$2
    par3=$3
    par4=$4
    ;;
5)  par1=$1
    par2=$2
    par3=$3
    par4=$4
    par5=$5
    ;;
6)  par1=$1
    par2=$2
    par3=$3
    par4=$4
    par5=$5
    par6=$6
    ;;
7)  par1=$1
    par2=$2
    par3=$3
    par4=$4
    par5=$5
    par6=$6
    par7=$7
    ;;
*) echo "Unsupported number of arguments" 
   echo "Usage : jobScript.sh par1 [par2] [par2] [par3] [par4] [par5] [par6] [par7]"
   exit
   ;;
esac

    echo ""               
    echo "--------------------------------"
    echo "SLURM_JOBID        : " $SLURM_JOBID
    echo "SLURM_ARRAY_JOB_ID : " $SLURM_ARRAY_JOB_ID
    echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
    echo "--------------------------------"
    echo ""

    arrayoffset=$par3 
    pathoutputlog=$par2 
    jobarrayFile=$par1

    ((myline=${SLURM_ARRAY_TASK_ID}+${arrayoffset}))
    # map back params for the job
    input=$(awk "NR==${myline}" $jobarrayFile)   # get all params for this job
    
    par1=$(echo $input | cut -d " " -f1)
    par2=$(echo $input | cut -d " " -f2)
    par3=$(echo $input | cut -d " " -f3)
    par4=$(echo $input | cut -d " " -f4)
    par5=$(echo $input | cut -d " " -f5)
    par6=$(echo $input | cut -d " " -f6)
    par7=$(echo $input | cut -d " " -f7)

    echo "par1 = ${par1}"
    echo "par2 = ${par2}"
    echo "par3 = ${par3}"
    echo "par4 = ${par4}"
    echo "par5 = ${par5}"
    echo "par6 = ${par6}"
    echo "par7 = ${par7}"

    format='+%Y/%m/%d-%H:%M:%S'

    host=$(hostname)

    date $format
    echo ""               
    echo "--------------------------------"
    echo "RUNNING ON HOST : " $host
    echo "WORKING DIR     : " $(pwd)
    echo "USER is         : " $USER
    echo "DISK USAGE /tmp :"
    df -h /tmp
    echo "--------------------------------" 
    
    echo ""               
    echo "--------------------------------"
    echo "par1 = ${par1}"
    echo "par2 = ${par2}"
    echo "par3 = ${par3}"
    echo "par4 = ${par4}"
    echo "par5 = ${par5}"
    echo "par6 = ${par6}"
    echo "par7 = ${par7}"
    echo "--------------------------------"
    echo ""

    echo ""
    echo "--------------------------------"
    echo " DEBUG INFO"                     
    echo "==> Kernel version information :"
    uname -a                               
    echo "==> C compiler that built the kernel:"
    cat /proc/version                           
    echo "==> load on this node:"               
    mpstat -P ALL                               
    echo "==> actual compiler is"               
    gcc -v                                      
    echo "--------------------------------"     
    echo ""               

###################################################################
###################################################################
#   EDIT THIS PART TO EXECUTE YOUR JOB!

#----------------------
# evironment 
  echo "==> running enironment script ${par1}"
. ${par1}
#----------------------

  echo "==> ldd ${par2}"
  ldd $par2

  echo "==> execute program "

  if [ $par7 == "single" ]
  then


    for file in $(echo $par3 | sed 's/,/ /g')
    do
        
        file_wo_path=$(basename $file | sed 's/.root//g')
        outfile_no_root=$(echo ${par4} | sed 's/.root//g')
        outfile=${outfile_no_root}_${file_wo_path}.root
        
        echo "==> $par2 $file ${outfile} $par5"
        time  $par2 $file ${outfile} $par5

        status=$?

        echo "------------------------------------"
        echo "OUTPUT:"

        outfile_wo_path=$(basename $outfile)

        if [ $status -ne 0 ]
        then
          echo "JOB: $JOB_ID CRASHED ON HOST: $host WITH OUTFILE ${outfile_wo_path}_${file_wo_path}"
          echo "JOB: $JOB_ID CRASHED ON HOST: $host WITH OUTFILE $outfile_wo_path_${file_wo_path}" > ${outdir}/crash/${outfile_wo_path}_${file_wo_path}_crash.txt
        fi
    done    
  #---------------------------------------------------------------------
  else

    echo "==> $par2 $par3 $par4 $par5"
    time  $par2 $par3 $par4 $par5

    status=$?

    echo "------------------------------------"
    echo "OUTPUT:"

    outfile=$par4
    inputfile=${par3}  # single file or comma separated list

    outfile_wo_path=$(basename $outfile)

    if [ $status -ne 0 ]
    then
      echo "JOB: $JOB_ID CRASHED ON HOST: $host WITH OUTFILE $outfile_wo_path"
      echo "JOB: $JOB_ID CRASHED ON HOST: $host WITH OUTFILE $outfile_wo_path" > ${outdir}/crash/${outfile_wo_path}_crash.txt
    fi

  fi

#   END EDIT YOUR EXECUT JOB!
###################################################################
###################################################################

    echo ""               
    echo "--------------------------------"
    echo "Job with params "
    echo "par1 = ${par1}"  
    echo "par2 = ${par2}"  
    echo "par3 = ${par3}"  
    echo "par4 = ${par4}"  
    echo "par5 = ${par5}"  
    echo "par6 = ${par6}"  
    echo "par7 = ${par7}"  
    echo "finsished!"      
    echo "--------------------------------"
    echo ""
    
    echo ""               
    echo "--------------------------------"
    echo "MONITOR ENVIRONMENT:"
    echo "SLURM_JOBID        : " $SLURM_JOBID
    echo "SLURM_ARRAY_JOB_ID : " $SLURM_ARRAY_JOB_ID
    echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
    echo "RUNNING ON HOST    : " $(hostname)
    echo "WORKING DIR        : " $(pwd)
    echo "USER is            : " $USER
    echo "DISK USAGE /tmp    :------------"
    df -h /tmp 
    echo "--------------------------------"
    
    date $format

    sleep 2

    outfile=$(basename $par4 | sed 's/.root//')

    mv ${pathoutputlog}/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${pathoutputlog}/${outfile}.log

# added propagation of the exit code from the root macro to here, so SLURM know if the job has crashed or not
    exit "$status"