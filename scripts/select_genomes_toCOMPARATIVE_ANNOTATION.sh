#!/bin/bash

# Check if the argument is provided
if [ $# -eq 0 ]; then
    echo "Please provide a search text as an argument."
    exit 1
fi

search_text=$1

output_dir="results/metagenome/BIN_ASSESSMENT/bins_quality_passedFinal"
mkdir -p "$output_dir"
head -n1 combined_metadata.csv > selected_metadata.csv
grep -w "$search_text" combined_metadata.csv >> selected_metadata.csv
while IFS=',' read -r genome rest; do
    # Skip the header row
    if [[ $genome == "Genome" ]]; then
        continue
    fi
    # Find the genome file and link it to the output directory
    genome_file=$(find results/metagenome/BIN_ASSESSMENT/bins_quality_passed_*/ -name "${genome}.fa")
    if [[ -n $genome_file ]]; then
        #ln -s "$genome_file" "$output_dir"
        cp "$genome_file" "$output_dir"
    fi
done < selected_metadata.csv

head -n1 selected_metadata.csv | tr "," "\n" | nl 
echo " " 
echo "##############################################################"
echo -e "$(ls ${PWD}/results/metagenome/BIN_ASSESSMENT/bins_quality_passedFinal/* | wc -l) genomes were selected " 
echo "Check selected_metadata.csv file and genomes in ${PWD}/results/metagenome/BIN_ASSESSMENT/bins_quality_passedFinal"

echo "### select metadata column and use it in COMPARATIVE_ANNOTATION"


echo " " 
echo "run below script"
echo "$ nextflow run COMPARATIVE_ANNOTATION.sif --metacol yourselectedcolumn"
echo "run above script" 
echo " " 
echo "You could choose one metadata column and run analyze specific interested genoems  "

 
