import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import os

def extract_gene_sequences(fasta_path, gff_path, selected_ids_path):
    # Get the base name of the input FASTA file without extension
    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    
    # Read the fasta file and create a dictionary of contigs
    contigs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    
    # Read the selected IDs and create a list of tuples (unique_id, full_description)
    selected_genes = []
    with open(selected_ids_path, "r") as file:
        for line in file:
            parts = line.strip().split('_')
            unique_id = f"{parts[1]}_{parts[2]}"
            selected_genes.append((unique_id, line.strip()))
    
    # Prepare to collect GFF data
    gff_data = {}
    with open(gff_path, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            contig, source, feature, start, end, score, strand, phase, attributes = parts
            if feature == "CDS":
                gene_id = attributes.split(';')[0].split('=')[1]
                gff_data[gene_id] = (contig, start, end, strand)
    
    # Extract sequences in the order of selected genes
    results = []
    for unique_id, full_desc in selected_genes:
        if unique_id in gff_data:
            contig, start, end, strand = gff_data[unique_id]
            # Extract sequence from contigs
            seq = contigs[contig].seq[int(start)-1:int(end)]
            if strand == '-':
                seq = seq.reverse_complement()
            # Create header with base name and full description
            header = f">{base_name}__{full_desc}"
            results.append(f"{header}\n{str(seq)}\n")
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Extract gene sequences based on selected gene IDs.")
    parser.add_argument("fasta_path", type=str, help="Path to the FASTA file containing contig sequences.")
    parser.add_argument("gff_path", type=str, help="Path to the GFF file containing gene annotations.")
    parser.add_argument("selected_ids_path", type=str, help="Path to the text file containing selected gene IDs.")

    args = parser.parse_args()

    # Run the function and capture the output
    output = extract_gene_sequences(args.fasta_path, args.gff_path, args.selected_ids_path)
    
    # Print the output
    for line in output:
        print(line)

if __name__ == "__main__":
    main()

