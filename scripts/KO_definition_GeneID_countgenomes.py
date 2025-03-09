import csv
from collections import defaultdict
import argparse

def load_ko_definitions(ko_list_file):
    ko_defs = {}
    with open(ko_list_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip header
        for row in reader:
            knum, definition = row[0], row[-1]
            # Replace commas with semicolons in the definition
            definition = definition.replace(',', ';')
            ko_defs[knum] = definition
    return ko_defs

def process_csv(input_file, output_file, ko_list_file):
    # Load KO definitions
    ko_defs = load_ko_definitions(ko_list_file)

    # Dictionary to store KO data
    ko_data = defaultdict(lambda: {'gene_ids': set(), 'genome_counts': defaultdict(int)})

    # Read input CSV
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            genome = row['Genome']
            gene_ids = row['GeneValue'].split()  # Split multiple Gene IDs
            ko = row['KO']
            if ko:  # Exclude rows with empty KO
                for gene_id in gene_ids:
                    ko_data[ko]['gene_ids'].add(gene_id)
                    ko_data[ko]['genome_counts'][genome] += 1

    # Write output CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        genomes = sorted(set(genome for ko_info in ko_data.values() for genome in ko_info['genome_counts']))
        writer.writerow(['KO', 'Definition', 'Gene IDs'] + genomes)

        # Write data rows
        for ko, info in ko_data.items():
            definition = ko_defs.get(ko, "")  # Get definition, or empty string if not found
            gene_ids = ';'.join(sorted(info['gene_ids']))
            counts = [info['genome_counts'].get(genome, 0) for genome in genomes]
            writer.writerow([ko, definition, gene_ids] + counts)

    print(f"Output written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Process KO gene data and add KEGG definitions.")
    parser.add_argument('-i', '--input', required=True, help="Input CSV file path")
    parser.add_argument('-o', '--output', required=True, help="Output CSV file path")
    parser.add_argument('-k', '--ko_list', required=True, help="KEGG KO list file path")
    
    args = parser.parse_args()

    process_csv(args.input, args.output, args.ko_list)

if __name__ == "__main__":
    main()
