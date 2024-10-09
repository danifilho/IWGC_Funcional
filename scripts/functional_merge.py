import argparse
import csv
import os
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Process annotation data.')
    parser.add_argument('--filelist', required=True, help='File containing list of input files')
    parser.add_argument('--output_dir', required=True, help='Directory to save output TSV files')
    parser.add_argument('--species_name', required=True, help='Species name used in file naming')
    args = parser.parse_args()

    # Initialize dictionaries
    ncbi_dict = defaultdict(list)
    ncbi_description_dict = {}
    signalp_dict = defaultdict(list)
    ipr_dict = defaultdict(set)
    go_dict = defaultdict(set)
    hmmer_dict = defaultdict(list)
    multiloc_dict = defaultdict(str)
    targetp_dict = defaultdict(str)
    other_mmseqs_data = defaultdict(dict)  # For other MMseqs databases

    # Collect all unique genes
    all_genes = set()

    # Read the list file and categorize input files
    mmseqs_files = []
    other_files = {
        'ncbi_mmseqs': '',
        'iprscan': '',
        'signalp': '',
        'targetp': '',
        'hmmer': '',
        'multiloc': '',
        'ncbidata': ''
    }

    with open(args.filelist, 'r') as filelist:
        for line in filelist:
            line = line.strip()
            basename = os.path.basename(line)
            if basename.endswith('_ncbi_best_hit.out'):
                other_files['ncbi_mmseqs'] = line
            elif basename.endswith('_best_hit.out'):
                mmseqs_files.append(line)
            elif basename.endswith('.iprscan.tsv'):
                other_files['iprscan'] = line
            elif basename.endswith('.gff3'):
                other_files['signalp'] = line
            elif basename.endswith('_summary.targetp2'):
                other_files['targetp'] = line
            elif basename.endswith('.domtblout'):
                other_files['hmmer'] = line
            elif basename.endswith('_ml2.txt'):
                other_files['multiloc'] = line
            elif basename.endswith('.prot'):
                other_files['ncbidata'] = line

    # Process NCBI descriptions
    print("Processing NCBI descriptions...")
    with open(other_files['ncbidata'], 'r') as ncbi_file:
        for line in ncbi_file:
            if line.startswith('>'):
                parts = line[1:].strip().split(' ', 1)
                accession = parts[0]
                description = parts[1] if len(parts) > 1 else "Description not found"
                ncbi_description_dict[accession] = description

    # Process NCBI MMseqs output
    print("Processing NCBI MMseqs data...")
    with open(other_files['ncbi_mmseqs'], 'r') as ncbi_mmseqs_file:
        for line in ncbi_mmseqs_file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                gene = fields[0]
                subject = fields[1]
                evalue = fields[10]
                bit_score = fields[11]
                ncbi_dict[gene] = [f"{subject};{evalue};{bit_score}"]
                all_genes.add(gene)

    # Multiloc
    print("Processing MultiLoc data...")
    with open(other_files['multiloc'], 'r') as multiloc:
        for line in multiloc:
            if not line.startswith('\t'):
                fields = line.strip().split()
                gene = fields[0]
                multiloc_dict[gene] = fields[1] + fields[2]
                all_genes.add(gene)

    # SignalP
    print("Processing SignalP data...")
    with open(other_files['signalp'], 'r') as signalp:
        for line in signalp:
            if not line.startswith('#'):
                fields = line.strip().split()
                gene = fields[0]
                signalp_dict[gene] = fields[4:6]
                all_genes.add(gene)

    # TargetP
    print("Processing TargetP data...")
    with open(other_files['targetp'], 'r') as targetpfile:
        for line in targetpfile:
            if not line.startswith('#'):
                fields = line.strip().split()
                gene = fields[0]
                targetp_dict[gene] = ';'.join(fields[1:])
                all_genes.add(gene)

    # Process other MMseqs outputs
    print("Processing other MMseqs data...")
    for file_path in mmseqs_files:
        if file_path == other_files['ncbi_mmseqs']:
            continue  # Skip the NCBI MMseqs file

        basename = os.path.basename(file_path)
        prefix = f"{args.species_name}_"
        suffix = "_best_hit.out"
        if basename.startswith(prefix) and basename.endswith(suffix):
            database_name = basename[len(prefix):-len(suffix)]
            if database_name == 'ncbi':
                continue  # Already processed NCBI MMseqs
        else:
            print(f"Error: MMseqs output file '{basename}' does not match expected pattern.")
            continue

        with open(file_path, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    fields = line.strip().split()
                    gene = fields[0]
                    all_genes.add(gene)
                    subject = fields[1]
                    evalue = fields[2]
                    bit_score = fields[3]
                    other_mmseqs_data[gene][database_name] = [subject, evalue, bit_score]


    # HMMER
    print("Processing HMMER data...")
    with open(other_files['hmmer'], 'r') as hmmerfile:
        for line in hmmerfile:
            if not line.startswith('#'):
                fields = line.strip().split()
                gene = fields[3]
                all_genes.add(gene)
                accession = fields[1]
                hmmer_dict[gene].append(accession)

    # IPRScan
    print("Processing IPRScan data...")
    with open(other_files['iprscan'], 'r') as iprscanfile:
        for line in iprscanfile:
            fields = line.strip().split('\t')
            gene = fields[0]
            all_genes.add(gene)
            ipr = fields[11] if fields[11] != '-' else ''
            go_terms = fields[13]
            if ipr:
                ipr_dict[gene].add(ipr)
            if go_terms:
                for go_term in go_terms.split('|'):
                    if 'GO:' in go_term:
                        go = go_term
                        go_dict[gene].add(go)

 # Build the output
    output_file = os.path.join(args.output_dir, f"functional_annotation.tsv")
    print("Writing results to:", output_file)

    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        # Build the header
        header = ["ID sequence", "NCBI DESCRIPTION", "NCBI Subject; E-value; Bit score"]

        # Get a sorted list of other MMseqs database names
        other_database_names = sorted({db for gene_data in other_mmseqs_data.values() for db in gene_data.keys()})

        # Add MMseqs columns dynamically for other databases
        for db_name in other_database_names:
            header.append(f"{db_name} Subject; E-value; Bit score")

        # Add the rest of the fixed columns
        header.extend([
            "SignalP Pos; Pr",
            "TargetP Prediction; noTP; SP; mTP; cTP; luTP; CS Position",
            "MULTILOC",
            "IPRSCAN GO",
            "IPRSCAN IPR",
            "Hmmer Pfam"
        ])
        writer.writerow(header)

        for gene in sorted(all_genes):
            row = [gene]
            # NCBI description and data
            ncbi_values = ncbi_dict.get(gene, ["N/A"])
            subject_id = ncbi_values[0].split(';')[0] if ncbi_values[0] != "N/A" else "N/A"
            description = ncbi_description_dict.get(subject_id, "Description not found")
            row.append(description)
            row.append(";".join(ncbi_values))

            # Other MMseqs data
            for db_name in other_database_names:
                db_values = other_mmseqs_data.get(gene, {}).get(db_name, ["N/A", "N/A", "N/A"])
                row.append(";".join(db_values))

            # Add fixed columns data
            row.extend([
                ";".join(signalp_dict.get(gene, ["N/A"])),
                targetp_dict.get(gene, "N/A"),
                multiloc_dict.get(gene, "N/A"),
                ";".join(sorted(go_dict.get(gene, ["N/A"]))),
                ";".join(sorted(ipr_dict.get(gene, ["N/A"]))),
                ";".join(hmmer_dict.get(gene, ["N/A"]))
            ])
            writer.writerow(row)

if __name__ == '__main__':
    main()
