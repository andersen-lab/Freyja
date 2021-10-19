bamfolder=$1
varfolder=$2
depthfolder=$3
# echo "${depthfolder}"
if [[ $# -eq 0 ]] ; then
    echo 'Need to supply input arguments.'
    exit 1
fi

for fn in $1/* ## bamfiles/* 
do
    # echo "${fn}"
    # echo "${fn%.*}" 
    # echo "${fn##*/}" 
    fn_out=${fn##*/}
    baseName=${fn##*/}
    baseName=${baseName%.*}
    echo $depthfolder$baseName
    samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f NC_045512_Hu-1.fasta ${fn} | tee >(cut -f1-4 > ${depthfolder}${fn_out%.*}.depth) | ivar variants -p $depthfolder$baseName -q 20 -r NC_045512_Hu-1.fasta
    # samtools mpileup -A -aa -d 600000 -Q 20 -q 0 -B -f NC_045512_Hu-1.fasta ${fn} | cut -f1-4 > ${depthfolder}${fn_out%.*}.depth #&
done


