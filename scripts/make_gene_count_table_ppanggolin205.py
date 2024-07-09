import os
import re
import argparse
from collections import defaultdict

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate gene count matrix from GFF files and gene_families.tsv')
parser.add_argument('input_dir', help='Input directory containing GFF folder and gene_families.tsv')
args = parser.parse_args()

# Extract genome ID from GFF files and create genome_id_linking file
genome_id_linking = {}
gff_dir = os.path.join(args.input_dir, 'gff')  # Set GFF directory path
for gff_file in os.listdir(gff_dir):
    with open(os.path.join(gff_dir, gff_file), "r") as f:
        for line in f:
            if "\tCDS\t" in line:
                match = re.search(r"ID=gnl\|Prokka\|(\w+)_(\d+)", line)
                #match = re.search(r"ID=gnl\|Prokka\|(\w+)_(\d+)", line)

                if match:
                    genome_id = match.group(1)
                    genome_id_linking[genome_id] = os.path.splitext(gff_file)[0]
                    break
with open(os.path.join(args.input_dir, "genome_id_linking.tsv"), "w") as f:
    for genome_id, gff_prefix in genome_id_linking.items():
        f.write(f"{gff_prefix}\t{genome_id}\n")

# Load gene families and associate with genome IDs
gene_families = defaultdict(lambda: defaultdict(int))
gene_families_file = os.path.join(args.input_dir, "gene_families.tsv")
with open(gene_families_file, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) >= 2:
            gene_family_id, gene_id = fields[0], fields[1]
            genome_id_prefix = re.search(r"gnl\|Prokka\|(.+?)_\d+", gene_id).group(1)
            genome_id = genome_id_linking.get(genome_id_prefix, genome_id_prefix)
            
            gene_families[gene_family_id][genome_id] += 1

# Output gene count matrix
gene_family_ids = sorted(gene_families.keys())
genome_ids = sorted(set(genome_id_linking.values()))

output_file = os.path.join(args.input_dir, "gene_count_matrix.tsv")
with open(output_file, "w") as f:
    f.write("Gene\t" + "\t".join(genome_ids) + "\n")
    for gene_family_id in gene_family_ids:
        counts = [str(gene_families[gene_family_id].get(genome_id, 0)) for genome_id in genome_ids]
        f.write(gene_family_id + "\t" + "\t".join(counts) + "\n")

