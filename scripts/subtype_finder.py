import os
import sys

def read_mapfile(mapfile):
    header_map = {}
    with open(mapfile, 'r') as file:
        for line in file:
            key, value = line.strip().split('|')
            header_map[key] = f'{key}|{value}'
    print(f"Mapfile read: {len(header_map)} entries found.")
    return header_map

def process_fasta_file(fasta_file, header_map):
    sequences = {}
    current_header = None
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line.strip().split(' ')[0][1:]  # Remove the '>'
                matched_key = None
                for key in header_map.keys():
                    if f'_{key}_' in f'_{header}_':
                        if matched_key is not None:
                            print(f"Warning: Multiple matches for header '{header}'. Skipping...")
                            matched_key = None
                            break
                        matched_key = key
                
                if matched_key:
                    new_header = header_map[matched_key]
                    current_header = new_header
                    if new_header not in sequences:
                        sequences[new_header] = []
                    sequences[new_header].append(f'>{new_header}\n')
                else:
                    current_header = None
            elif current_header:
                sequences[current_header].append(line)
    print(f"Processed {fasta_file}: {len(sequences)} matching headers found.")
    return sequences

def write_sequences_to_files(sequences, result_dir):
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    for header, seq_lines in sequences.items():
        subtype = header.split('|')[1]
        output_file = os.path.join(result_dir, f'{subtype}.fna')
        with open(output_file, 'a') as file:
            file.writelines(seq_lines)
    print(f"Written sequences to {result_dir}")

def process_directory(fasta_dir, mapfile, result_dir):
    header_map = read_mapfile(mapfile)
    for filename in os.listdir(fasta_dir):
        if filename.endswith('_GHoINT.ffn'):
            fasta_file = os.path.join(fasta_dir, filename)
            sequences = process_fasta_file(fasta_file, header_map)
            write_sequences_to_files(sequences, result_dir)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <fasta_directory> <mapfile> <result_directory>")
        sys.exit(1)

    fasta_directory = sys.argv[1]
    mapfile = sys.argv[2]
    result_directory = sys.argv[3]

    print(f"Processing fasta directory: {fasta_directory}")
    print(f"Using mapfile: {mapfile}")
    print(f"Results will be saved in: {result_directory}")

    process_directory(fasta_directory, mapfile, result_directory)
    print("Processing completed.")

