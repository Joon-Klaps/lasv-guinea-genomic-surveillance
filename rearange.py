#! /usr/bin/env python3

def main ():
    """
    Small python script that rearanges the fasta file sequences with a given order.
    Usage: python3 rearange.py <fasta_file> <order_file> <output_file>
    """
    import sys

    fasta_file = sys.argv[1]
    order_file = sys.argv[2]
    output_file = sys.argv[3]

    # Read the order of sequence IDs
    with open(order_file, 'r') as f:
        order = [line.strip() for line in f if line.strip()]

    # Read the fasta file and store sequences in a dictionary
    sequences = {}
    with open(fasta_file, 'r') as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    sequences[seq_id] = ''.join(seq_lines)
                seq_id = line[1:]  # Remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id is not None:
            sequences[seq_id] = ''.join(seq_lines)

    # Write the sequences in the specified order to the output file
    with open(output_file, 'w') as f:
        for seq_id in order:
            if seq_id in sequences:
                f.write(f'>{seq_id}\n')
                f.write(f'{sequences[seq_id]}\n')
                sequences.pop(seq_id)
        for seq_id in sequences:
                f.write(f'>{seq_id}\n')
                f.write(f'{sequences[seq_id]}\n')

if __name__ == "__main__":
    main()