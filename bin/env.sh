bam=$1
bed=$2
ref=$3
output=$4
export PYTHON=!{params.python}
OUTPUT_FOLDER=$(realpath ${output})
export OUTPUT_FOLDER
echo "$OUTPUT_FOLDER"
BAM=$(realpath ${bam})
export BAM
BED=$(realpath ${bed})
export BED
REF=$(realpath ${ref})
export REF
LOG_PATH="${OUTPUT_FOLDER}/log"
export LOG_PATH