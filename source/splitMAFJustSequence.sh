if [[ $# -ne 5 ]]; then
    echo 'Usage: ./splitMAFJustSequence.sh <MAF sequence file> <split size> <output folder> <output prefix> <W parameter>'
    exit 0
fi

MAFfile=$1
splitSize=$2
odir=$3
oprefix=$4
W=$5

mkdir -p ${odir}

split -l ${splitSize} -d <(zcat ${MAFfile}) ${odir}/${oprefix}_maf_sequence

numFiles=`wc -l <(ls -l ${odir}/${oprefix}_maf_sequence) |cut -f1 -d ' '`
firstLine=`head -n 1 ${odir}/${oprefix}_maf_sequence00`
for ((i=0; i<$numFiles; i++))
do
    cur="${odir}/${oprefix}_maf_sequence$(printf "%02d" $i)"
    next="${odir}/${oprefix}_maf_sequence$(printf "%02d" $((i+1)))"

    if [ -f "$next" ]; then
        tail -n ${W} ${cur} > temp_${oprefix}.txt # put last 10 in a temp file
        head -n ${W} ${next} >> ${cur} # add first 10 lines of next file to end of current file
        echo ${firstLine} > temp_next_${oprefix}.txt
        cat temp_${oprefix}.txt ${next} >> temp_next_${oprefix}.txt # add last 10 lines of current file to beginning of next 
        mv temp_next_${oprefix}.txt ${next}
    fi

    gzip ${cur}
done
