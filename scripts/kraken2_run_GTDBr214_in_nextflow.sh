#!/bin/bash 
set -e 
set -u 

database_dir=${1}
accession=${2}
cpus=${3}
cpus=$((cpus + 0))
read1=${4}
read2=${5}

time kraken2 --db ${database_dir} --output ${accession}_kraken2.out \
--report ${accession}_kraken2.report --confidence 0.01 --gzip-compressed --paired \
--threads ${cpus} ${read1} ${read2}

# --memory-mapping 
