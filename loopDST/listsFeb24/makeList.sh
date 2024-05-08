

if [ $# -eq 1 ]
then
    indir=$1

    dayold=""

    for dir in $(find $indir -maxdepth 1 -type d)
    do
        
        day=$(basename $dir)
        
        if [ "$dir" != "$indir" ]
        then
            outfile=day_${day}.list 
        
            if [ -e $outfile ]
            then
                rm -f $outfile
            fi
            
            echo "NEW FILE $day"
            
            for i in {1..8}
            do
                find $dir/0${i}/root -maxdepth 1 -type f -name "be*_dst_feb24*.root"  >> $outfile
            done
            nfiles=$(cat $outfile  | wc -l)
            echo "FILE ${outfile} contains ${nfiles} files"
        fi

    done

else
  echo  "input dir needed!"
fi  
