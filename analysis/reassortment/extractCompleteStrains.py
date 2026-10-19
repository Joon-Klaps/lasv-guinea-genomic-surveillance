# !/usr/bin/env python3

import pandas as pd
import argparse
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Dict, List, Optional
import re

def main():
    parser = argparse.ArgumentParser(description='Extract complete strains from a metadata file.')
    parser.add_argument('--xlsx',
                        help='Path to the metadata file (Excel format).',
                        default="../../data/curated-metadata/LASV_GENBANK_20250725_curated.raw.xlsx")
    parser.add_argument('--xml_files',
                        nargs='+',
                        help='Paths to the beast XML files.',
                        default=["../../analysis/phylogeographic/runs/L/LASV_GENBANK_20250725_curated.L.noNIG-CIV.beast2.xml",
                                "../../analysis/phylogeographic/runs/S/LASV_GENBANK_20250725_curated.S.noNIG-CIV.beast2.xml"])
    parser.add_argument('--prefix',
                        help='Prefix for the output alignment files.',
                        default="../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed")
    args = parser.parse_args()

    sequences = {}
    for xml_file in args.xml_files:
        sequences.update(xml_to_seq(xml_file))

    # Read the metadata file
    df = pd.read_excel(args.xlsx)

    # Filter for complete strains
    complete_strains = identify_complete_strains(df, sequences)

    # Save the complete strains to the output file
    for segment, seq_records in complete_strains.items():
        output_file = f"{args.prefix}.{segment}.msa"
        SeqIO.write(seq_records, output_file, "fasta")
        print(f"Saved {len(seq_records)} complete strains for segment {segment} to {output_file}")

def xml_to_seq(xml_file: str) -> Dict[str, SeqRecord]:
    """
    Extract sequences from a BEAST XML file.

    Region to look for:
        <alignment id="alignment" dataType="nucleotide">
		<sequence>
			<taxon idref="MT119509|Dam-67|Mastomys_natalensis|GIN|2014-03-26"/>
            ACTG...
		</sequence>

    Args:
        xml_file (str): Path to the BEAST XML file.

    Returns:
        Dict[str, SeqIO.SeqRecord]: Sequence records extracted from the XML file.
    """
    sequences: Dict[str, SeqRecord] = {}

    tree = ET.parse(xml_file)
    root = tree.getroot()

    for seq_el in root.findall(".//alignment//sequence"):
        taxon_el = seq_el.find("taxon")
        if taxon_el is None:
            continue

        taxon_id = taxon_el.attrib.get("idref") or taxon_el.attrib.get("id")
        if not taxon_id:
            continue

        seq_text_parts = [text.strip() for text in seq_el.itertext() if text and text.strip()]
        if not seq_text_parts:
            continue

        sequence = "".join(seq_text_parts).replace(" ", "").replace("\n", "").upper()
        sequences[taxon_id] = SeqRecord(seq=Seq(sequence), id=taxon_id, description="")

    return sequences




def identify_complete_strains(df, sequences: Dict[str, SeqRecord]) -> Dict[str, List[SeqRecord]]:
    """
    Identify complete strains from the metadata dataframe.
    """
    # Normalize Strain
    df = df.copy()
    df['Strain_normalized'] = df['Strain'].astype(str).str.lower().str.replace(r'(.*-\d)-[12]$', r'\1', regex=True)

    # 1. Create Lookup Maps
    # Using dropna and set_index to create fast lookup dictionaries.
    # Prioritizing the first occurrence for duplicates as per original logic.
    df_align = df.dropna(subset=['alignment_name']).set_index('alignment_name', drop=False)
    df_align = df_align[~df_align.index.duplicated(keep='first')]

    df_acc = df.dropna(subset=['Accession']).set_index('Accession', drop=False)
    df_acc = df_acc[~df_acc.index.duplicated(keep='first')]

    matched_rows = []
    unmatched_sequences = []

    # 2. Match Sequences to Metadata
    for seq_id in sequences:
        row = None
        # Try alignment_name
        if seq_id in df_align.index:
            row = df_align.loc[seq_id]
        else:
            # Try Accession
            accession = seq_id.split("|")[0] if "|" in seq_id else seq_id
            if accession in df_acc.index:
                row = df_acc.loc[accession]

        if row is not None:
             r = row.to_dict()
             r['__seq_id'] = seq_id
             matched_rows.append(r)
        else:
            segment = 'Unknown'
            if '.L.' in seq_id or seq_id.endswith('_L'):
                segment = 'L'
            elif '.S.' in seq_id or seq_id.endswith('_S'):
                segment = 'S'

            unmatched_sequences.append({
                'alignment_name': seq_id,
                'extracted_accession': seq_id.split("|")[0] if "|" in seq_id else seq_id,
                'segment': segment,
                'reason': 'Not found in metadata'
            })

    # 3. Handle Unmatched
    if unmatched_sequences:
        df_unmatched = pd.DataFrame(unmatched_sequences)
        tmp_file = "sequences_not_in_dataframe.tsv"
        df_unmatched.to_csv(tmp_file, sep='\t', index=False)
        print(f"Saved {len(unmatched_sequences)} sequences not in dataframe to {tmp_file}")
        print(f"  L segments: {len(df_unmatched[df_unmatched['segment'] == 'L'])}")
        print(f"  S segments: {len(df_unmatched[df_unmatched['segment'] == 'S'])}")
        print(f"  Unknown: {len(df_unmatched[df_unmatched['segment'] == 'Unknown'])}")

    if not matched_rows:
        print("Warning: No sequences matched to dataframe!")
        return {"L": [], "S": []}

    matched = pd.DataFrame(matched_rows)
    complete = {"L": [], "S": []}
    removed_log = []

    # 4. Group by Strain and Validate
    for strain_norm, group in matched.groupby('Strain_normalized'):
        segments = set(group['Segment'].astype(str))

        if segments == {'L', 'S'}:
            # Has both segments
            for segment in ['L', 'S']:
                seg_rows = group[group['Segment'] == segment]

                # Take first, log others as duplicates
                target = seg_rows.iloc[0]
                if len(seg_rows) > 1:
                     print(f"Warning: Multiple {segment} segments for strain {strain_norm}, using first")
                     for _, dup in seg_rows.iloc[1:].iterrows():
                         removed_log.append({
                             **dup.to_dict(),
                             'strain': strain_norm,
                             'reason': 'Duplicate segment for strain',
                             'segments_found': 'L, S'
                         })

                seq_id = target['__seq_id']
                if seq_id in sequences:
                    rec = sequences[seq_id]
                    # Create new record with Strain as ID
                    new_rec = SeqRecord(rec.seq, id=strain_norm, description=rec.id)
                    complete[segment].append(new_rec)
                else:
                    # Should not happen given how we built matched
                    removed_log.append({
                        **target.to_dict(),
                        'strain': strain_norm,
                        'reason': 'Sequence missing in map (unexpected)',
                        'segments_found': 'L, S'
                    })
        else:
             # Incomplete strain
             reason = f"Incomplete: has segments {segments}"
             for _, row in group.iterrows():
                 if str(row['Segment']) in ['L', 'S']:
                     removed_log.append({
                         **row.to_dict(),
                         'strain': strain_norm,
                         'reason': reason,
                         'segments_found': ', '.join(sorted(segments))
                     })

    # 5. Save Removed Log
    if removed_log:
        df_rem = pd.DataFrame(removed_log)
        for segment in ['L', 'S']:
             rem_subset = df_rem[df_rem['Segment'] == segment]
             if not rem_subset.empty:
                tmp_file = f"removed_{segment}_sequences.tsv"
                rem_subset.to_csv(tmp_file, sep='\t', index=False)
                print(f"Saved {len(rem_subset)} removed {segment} sequences to {tmp_file}")

    return complete

if __name__ == "__main__":
    main()
