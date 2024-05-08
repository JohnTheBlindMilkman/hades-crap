#!/bin/bash

currentDir=$(pwd | xargs -i basename {})
currentDir=../$currentDir
prefix=6sec

#DAY=096
#INDIR=/lustre/nyx/hades/raw/apr12/${DAY}/
OUTDIR=${currentDir}/lists/${prefix}

#ls -A ${INDIR} > ${OUTFILE}

#suggested by jochen
selection=/lustre/hades/dst/apr12/gen8/sector_selection/FileListHadron.list
dstdir=/lustre/hades/dst/apr12/gen8
rawdir=/lustre/hades/raw/apr12

#seperated into 6-sector and 5-sector files, change manually!
#for file in $(cat $selection | grep "1 1 0 1 1 1" | cut -d " " -f 1)
for file in $(cat $selection | grep "1 1 1 1 1 1" | cut -d " " -f 1)
do
    day=${file:4:3}
    #find ${dstdir}/${day}/root/${file}*
    ls ${rawdir}/${day}/${file}* >> ${OUTDIR}/day_${day}.list
    
done
