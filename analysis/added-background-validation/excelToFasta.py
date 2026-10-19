#!/usr/bin/env python3
"""
Modular LASV sequence processing tool for converting Excel data to FASTA files,
generating YAML configurations, and cleaning alignment and tree files.
"""

import pandas as pd
import argparse
import yaml
from pathlib import Path
from typing import Dict, Set, Optional
import logging
from Bio import SeqIO, Phylo

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LASVDataProcessor:
    """Main class for processing LASV sequence data from Excel files."""

    def __init__(self, input_file: str, output_prefix: str):
        """
        Initialize the processor with input and output paths.

        Args:
            input_file: Path to input Excel file
            output_prefix: Prefix for output files
        """
        self.input_file = Path(input_file)
        self.output_prefix = output_prefix
        self.df: Optional[pd.DataFrame] = None
        self.filtered_df: Optional[pd.DataFrame] = None

    def load_data(self) -> pd.DataFrame:
        """Load data from Excel file."""
        logger.info("Loading data from %s", self.input_file)
        self.df = pd.read_excel(self.input_file)
        logger.info("Loaded %d rows from Excel file", len(self.df))
        return self.df

    def filter_data(self) -> pd.DataFrame:
        """Filter data based on alignment conditions."""
        if self.df is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        logger.debug("Starting data filtering process")
        logger.debug("Original dataset size: %d rows", len(self.df))

        # Check what columns exist
        logger.debug("Available columns: %s", list(self.df.columns))

        kept_condition = self.df['kept_in_alignment'].fillna(True).astype(bool)
        alignment_condition = self.df['last_alignment'].fillna(0) != 11

        logger.debug("Sequences with kept_in_alignment=True: %d", kept_condition.sum())
        logger.debug("Sequences with last_alignment!=11: %d", alignment_condition.sum())

        filtered_result = self.df[kept_condition & alignment_condition]
        # Ensure we have a DataFrame
        if isinstance(filtered_result, pd.Series):
            # This shouldn't happen with the current logic, but just in case
            raise ValueError("Unexpected Series result from DataFrame filtering")

        self.filtered_df = filtered_result

        logger.info("After filtering: %d rows", len(self.filtered_df))
        logger.debug("Filtered out %d rows", len(self.df) - len(self.filtered_df))

        # Add segment breakdown
        if 'Segment' in self.filtered_df.columns:
            segment_counts = self.filtered_df['Segment'].value_counts()
            for segment, count in segment_counts.items():
                logger.debug("Segment %s: %d sequences", segment, count)

        return self.filtered_df

    def get_sequences_to_remove(self) -> Set[str]:
        """Get set of accession numbers that should be removed from starting alignment."""
        if self.df is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        to_remove_condition = self.df['to_remove_from_starting_alignment'].fillna(False).astype(bool)
        to_remove = self.df[to_remove_condition]

        # Get original accessions and ensure they are strings
        original_accessions = {str(acc) for acc in to_remove['Accession']}

        # Create variants with "|" replaced by "/"
        # This handles cases where alignment tools have modified the accession format
        all_accessions = set(original_accessions)
        modified_accessions1 = [accession.replace("|", "/") for accession in original_accessions if "|" in accession]
        modified_accessions2 = [accession.replace("|", "_") for accession in original_accessions if "|" in accession]
        all_accessions.update(modified_accessions1)
        all_accessions.update(modified_accessions2)

        logger.info("Found %d sequences to remove from starting alignment", len(to_remove))
        logger.info("Generated %d total accession variants (including originals)", len(all_accessions))
        logger.debug("Accessions to remove: %s", all_accessions)
        return all_accessions

    @staticmethod
    def get_existing_sequence_ids(alignment_file: str) -> Set[str]:
        """
        Extract sequence IDs from an existing alignment file.

        Args:
            alignment_file: Path to the alignment file

        Returns:
            Set of sequence IDs found in the alignment
        """
        if not alignment_file or not Path(alignment_file).exists():
            logger.warning("Alignment file %s does not exist", alignment_file)
            return set()

        logger.debug("Extracting sequence IDs from %s", alignment_file)
        sequence_ids = set()

        try:
            with open(alignment_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        seq_id = line[1:]  # Remove '>'
                        sequence_ids.add(seq_id)

            logger.info("Extracted %d sequence IDs from %s", len(sequence_ids), alignment_file)
            return sequence_ids

        except IOError as e:
            logger.error("Error reading alignment file %s: %s", alignment_file, e)
            return set()
class YAMLConfigGenerator:
    """Handles generation of YAML configuration files."""

    def __init__(self):
        self.config = {
            'np': {'to_remove': [], 'to_remove_chunk': {'pattern': 'nnnnnnnnn', 'leading': [], 'trailing': []}},
            'gpc': {'to_remove': [], 'to_remove_chunk': {'pattern': 'nnnnnnnnn', 'leading': [], 'trailing': []}},
            'z': {'to_remove': [], 'to_remove_chunk': {'pattern': 'nnnnnnnnn', 'leading': [], 'trailing': []}},
            'l': {'to_remove': [], 'to_remove_chunk': {'pattern': 'nnnnnnnnn', 'leading': [], 'trailing': []}}
        }

    def process_sequences(self, df: pd.DataFrame) -> Dict:
        """Process sequences and populate YAML config."""
        logger.debug("Starting YAML config generation")
        logger.debug("Processing %d sequences for YAML config", len(df))

        gene_counts = {'GPC': 0, 'NP': 0, 'Z': 0, 'L': 0, 'GPC-NP': 0, 'Z-L': 0, 'other': 0}

        for _, row in df.iterrows():
            accession = str(row['Accession'])
            name = str(row['alignment_name'])
            genes = str(row['genes'])

            logger.debug("Processing sequence: %s (accession: %s, genes: %s)", name, accession, genes)

            if genes == 'GPC':
                self.config['np']['to_remove'].append(name)
                gene_counts['GPC'] += 1
            elif genes == 'NP':
                self.config['gpc']['to_remove'].append(name)
                gene_counts['NP'] += 1
            elif genes == 'Z':
                self.config['l']['to_remove'].append(name)
                gene_counts['Z'] += 1
            elif genes == 'L':
                self.config['z']['to_remove'].append(name)
                gene_counts['L'] += 1
            elif '-' in accession and genes == 'GPC-NP':
                self.config['gpc']['to_remove_chunk']['trailing'].append(name)
                self.config['np']['to_remove_chunk']['leading'].append(name)
                gene_counts['GPC-NP'] += 1
                logger.debug("Added %s to GPC trailing and NP leading removal chunks", name)
            elif '-' in accession and genes == 'Z-L':
                self.config['z']['to_remove_chunk']['trailing'].append(name)
                self.config['l']['to_remove_chunk']['leading'].append(name)
                gene_counts['Z-L'] += 1
                logger.debug("Added %s to Z trailing and L leading removal chunks", name)
            else:
                gene_counts['other'] += 1
                logger.debug("Unknown gene type: %s for sequence %s", genes, name)

        # Log summary of gene types processed
        for gene_type, count in gene_counts.items():
            if count > 0:
                logger.info("Processed %d sequences with gene type: %s", count, gene_type)

        return self.config

    def save_config(self, output_path: str) -> None:
        """Save YAML configuration to file."""
        with open(output_path, 'w', encoding='utf-8') as f:
            yaml.dump(self.config, f, default_flow_style=False)
        logger.info("YAML configuration saved to %s", output_path)


class FastaGenerator:
    """Handles FASTA file generation."""

    @staticmethod
    def save_sequences(df: pd.DataFrame, output_prefix: str,
                      alignment_s: Optional[str] = None,
                      alignment_l: Optional[str] = None) -> None:
        """
        Save sequences to FASTA files by segment, excluding sequences already in base alignments.

        Args:
            df: DataFrame containing sequences to save
            output_prefix: Prefix for output files
            alignment_s: Path to S segment base alignment (sequences to exclude)
            alignment_l: Path to L segment base alignment (sequences to exclude)
        """
        logger.debug("Starting to save sequences with output prefix: %s", output_prefix)

        # Get existing sequence IDs from base alignments
        existing_s_ids = set()
        existing_l_ids = set()

        if alignment_s:
            existing_s_ids = LASVDataProcessor.get_existing_sequence_ids(alignment_s)
            logger.info("Found %d existing sequences in S alignment", len(existing_s_ids))

        if alignment_l:
            existing_l_ids = LASVDataProcessor.get_existing_sequence_ids(alignment_l)
            logger.info("Found %d existing sequences in L alignment", len(existing_l_ids))

        for segment in ['L', 'S']:
            seg_data = df[df['Segment'] == segment]
            logger.debug("Processing %d sequences for segment %s", len(seg_data), segment)

            # Determine which existing IDs to exclude for this segment
            existing_ids = existing_s_ids if segment == 'S' else existing_l_ids

            # Filter out sequences that already exist in base alignment
            filtered_sequences = []
            excluded_count = 0

            for _, row in seg_data.iterrows():
                alignment_name = str(row['alignment_name'])
                accession = str(row['Accession'])

                # Check if this sequence is already in the base alignment
                is_duplicate = False
                for existing_id in existing_ids:
                    if accession in existing_id or alignment_name in existing_id:
                        is_duplicate = True
                        excluded_count += 1
                        logger.debug("Excluding sequence %s (accession: %s) - already exists in base alignment as %s",
                                   alignment_name, accession, existing_id)
                        break

                if not is_duplicate:
                    filtered_sequences.append(row)

            # Write filtered sequences to file
            output_file = f"{output_prefix}_{segment}.fasta"

            with open(output_file, 'w', encoding='utf-8') as f:
                for row in filtered_sequences:
                    f.write(f">{row['alignment_name']}\n{row['sequence']}\n")

            logger.info("Saved %d NEW sequences to %s (excluded %d duplicates from base alignment)",
                       len(filtered_sequences), output_file, excluded_count)

            if excluded_count > 0:
                logger.info("Segment %s: Excluded %d sequences already present in base alignment",
                           segment, excluded_count)
class AlignmentCleaner:
    """Handles cleaning of alignment files by removing specified sequences."""

    @staticmethod
    def clean_alignment(input_file: str, output_file: str, accessions_to_remove: Set[str]) -> None:
        """
        Remove sequences from alignment file based on exact accession matching.
        Avoids removing sequences where the accession is part of a range (e.g., MT119464-MT119558).

        Args:
            input_file: Path to input alignment file
            output_file: Path to output cleaned alignment file
            accessions_to_remove: Set of accession numbers to remove
        """
        if not Path(input_file).exists():
            logger.warning("Alignment file %s does not exist, skipping", input_file)
            return

        logger.info("Cleaning alignment file: %s", input_file)

        with open(input_file, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8') as outfile:
            current_seq_id = None
            current_seq_lines = []
            should_keep = True
            removed_count = 0
            kept_count = 0

            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous sequence if any
                    if current_seq_id is not None:
                        if should_keep:
                            outfile.write(f">{current_seq_id}\n")
                            outfile.write('\n'.join(current_seq_lines) + '\n')
                            kept_count += 1
                        else:
                            removed_count += 1

                    # Start new sequence
                    current_seq_id = line[1:]  # Remove '>'
                    current_seq_lines = []

                    # Check if we should keep this sequence
                    should_keep = True
                    for accession in accessions_to_remove:
                        if accession in current_seq_id:
                            # Accomodate for the collapsed sequence ID MK00001 in MK00001-MK00002 should not be removed
                            has_dash_before = '-' + accession in current_seq_id
                            has_dash_after = accession + '-' in current_seq_id
                            # Only remove if accession appears standalone (not as part of a range)
                            if not has_dash_before and not has_dash_after:
                                should_keep = False
                                logger.debug("Removing sequence %s (matches %s)", current_seq_id, accession)
                                break
                else:
                    current_seq_lines.append(line)

            # Process last sequence
            if current_seq_id is not None:
                if should_keep:
                    outfile.write(f">{current_seq_id}\n")
                    outfile.write('\n'.join(current_seq_lines) + '\n')
                    kept_count += 1
                else:
                    removed_count += 1

        logger.info("Cleaned alignment saved to %s", output_file)
        logger.info("Removed %d sequences, kept %d", removed_count, kept_count)
class TreeCleaner:
    """Handles cleaning of phylogenetic tree files by removing specified taxa."""

    @staticmethod
    def clean_tree(input_file: str, output_file: str, accessions_to_remove: Set[str],
                   alignment_file: Optional[str] = None) -> None:
        """
        Remove taxa from tree file based on exact accession matching and alignment consistency.
        Avoids removing taxa where the accession is part of a range (e.g., MT119464-MT119558).

        Args:
            input_file: Path to input tree file
            output_file: Path to output cleaned tree file
            accessions_to_remove: Set of accession numbers to remove
            alignment_file: Path to corresponding alignment file (if provided, removes sequences
                          that are in tree but missing from alignment)
        """

        if not Path(input_file).exists():
            logger.warning("Tree file %s does not exist, skipping", input_file)
            return

        logger.info("Cleaning tree file: %s", input_file)

        try:
            # Read tree
            tree = Phylo.read(input_file, "newick")
            original_terminals = len(tree.get_terminals())

            # Get sequences present in alignment if alignment file is provided
            alignment_sequences = set()
            if alignment_file and Path(alignment_file).exists():
                alignment_sequences = LASVDataProcessor.get_existing_sequence_ids(alignment_file)
                logger.info("Found %d sequences in alignment file %s", len(alignment_sequences), alignment_file)
            elif alignment_file:
                logger.warning("Alignment file %s does not exist, skipping alignment consistency check", alignment_file)

            # Get all terminal names
            all_terminals = tree.get_terminals()
            terminal_names_to_remove = set()

            # Find terminals that match removal list using smart matching logic
            for terminal in all_terminals:
                terminal_name = terminal.name or ""
                should_remove = False

                # Check against accessions_to_remove with smart matching
                for accession in accessions_to_remove:
                    if accession in terminal_name:
                        has_dash_before = '-' + accession in terminal_name
                        has_dash_after = accession + '-' in terminal_name

                        # Only remove if accession appears standalone (not as part of a range)
                        if not has_dash_before and not has_dash_after:
                            should_remove = True
                            terminal_names_to_remove.add(terminal_name)
                            logger.debug("Marking terminal %s for removal (matches %s)", terminal_name, accession)
                            break

            # Also remove terminals that are in tree but missing from alignment
            if alignment_file and alignment_sequences:
                all_terminal_names = {terminal.name or "" for terminal in all_terminals}
                missing_from_alignment = all_terminal_names - alignment_sequences

                if missing_from_alignment:
                    logger.info("Found %d terminals in tree that are missing from alignment", len(missing_from_alignment))
                    for missing_name in missing_from_alignment:
                        terminal_names_to_remove.add(missing_name)
                        logger.debug("Marking terminal %s for removal (missing from alignment)", missing_name)

            # Check for accessions that weren't found in the tree
            tree_terminal_names = {terminal.name or "" for terminal in all_terminals}
            unmatched_accessions = []
            for accession in accessions_to_remove:
                found_match = False
                for terminal_name in tree_terminal_names:
                    if accession in terminal_name:
                        found_match = True
                        break
                if not found_match:
                    unmatched_accessions.append(accession)

            if unmatched_accessions:
                logger.warning("The following accessions were not found in the tree: %s",
                               ", ".join(unmatched_accessions))

            # Remove terminals
            removed_count = 0
            if terminal_names_to_remove:
                logger.info("Removing %d terminals from tree", len(terminal_names_to_remove))
                for terminal_name in terminal_names_to_remove:
                    try:
                        # Find the terminal object by name and remove it
                        for terminal in tree.get_terminals():
                            if terminal.name == terminal_name:
                                tree.prune(terminal)
                                removed_count += 1
                                logger.debug("Removed terminal: %s", terminal_name)
                                break
                    except Exception as e:
                        logger.warning("Failed to remove terminal %s: %s", terminal_name, e)
            else:
                logger.info("No terminals matched for removal, tree remains unchanged")

            # Write cleaned tree
            Phylo.write(tree, output_file, "newick")
            final_terminals = len(tree.get_terminals())

            logger.info("Cleaned tree saved to %s", output_file)
            logger.info("Removed %d terminals, kept %d out of %d", removed_count, final_terminals, original_terminals)

        except Exception as e:
            logger.error("Error processing tree file %s: %s", input_file, e)


def get_cleaned_output_path(input_path_str: str) -> str:
    """
    Generate output path for cleaned files.
    If input is in 'archive' directory, save to parent directory.
    Otherwise, save in same directory as input.

    Args:
        input_path_str: Input file path as string

    Returns:
        Output file path with .cleaned suffix before extension
    """
    input_path = Path(input_path_str)

    if "archive" in input_path.parts:
        parent_dir = input_path.parent.parent
        output_path = parent_dir / f"{input_path.stem}.cleaned{input_path.suffix}"
    else:
        output_path = input_path.parent / f"{input_path.stem}.cleaned{input_path.suffix}"

    return str(output_path)


def main():
    parser = argparse.ArgumentParser(description='Convert Excel to FASTA files and clean alignments/trees')
    parser.add_argument('--input', default="../../data/curated-metadata/LASV_GENBANK_20250725_curated.manual.xlsx",
                        help='Input Excel file path')
    parser.add_argument('--output', default="../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.new-guinea-sequences",
                        help='Output file prefix (will create _L.fasta and _S.fasta)')
    parser.add_argument('--alignment_s',default="../../data/alignments/starting-alignment/archive/LASV_GENBANK_20250725_curated.background_S.msa",
                        help='Path to S segment alignment file to clean')
    parser.add_argument('--alignment_l', default="../../data/alignments/starting-alignment/archive/LASV_GENBANK_20250725_curated.background_L.msa",
                        help='Path to L segment alignment file to clean')
    parser.add_argument('--tree_s',default="../../data/trees/archive/LASV_GENBANK_20250725_curated.added_background_S_concat.fa.treefile",
                        help='Path to S segment tree file to clean')
    parser.add_argument('--tree_l', default="../../data/trees/archive/LASV_GENBANK_20250725_curated.added_background_L_concat.fa.treefile",
                        help='Path to L segment tree file to clean')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug logging for detailed output')

    args = parser.parse_args()

    # Set logging level based on debug flag
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info("Debug logging enabled")

    logger.info("Starting LASV data processing pipeline")
    logger.debug("Input file: %s", args.input)
    logger.debug("Output prefix: %s", args.output)
    logger.debug("S alignment: %s", args.alignment_s)
    logger.debug("L alignment: %s", args.alignment_l)

    # Initialize processor
    processor = LASVDataProcessor(args.input, args.output)

    # Load and filter data
    processor.load_data()
    filtered_data = processor.filter_data()

    # Generate YAML config
    yaml_generator = YAMLConfigGenerator()
    yaml_generator.process_sequences(filtered_data)
    yaml_generator.save_config(f"{args.output}_filter.yml")

    # Generate FASTA files (excluding sequences already in base alignments)
    FastaGenerator.save_sequences(filtered_data, args.output,
                                 args.alignment_s, args.alignment_l)

    # Get sequences to remove from starting alignment
    sequences_to_remove = processor.get_sequences_to_remove()

    if sequences_to_remove:
        logger.info("Processing removal of %d sequences from alignments and trees", len(sequences_to_remove))
        logger.debug("Sequences to remove: %s", list(sequences_to_remove))

        # Clean alignment files
        if args.alignment_s:
            logger.debug("Cleaning S segment alignment: %s", args.alignment_s)
            AlignmentCleaner.clean_alignment(
                args.alignment_s,
                get_cleaned_output_path(args.alignment_s),
                sequences_to_remove
            )

        if args.alignment_l:
            logger.debug("Cleaning L segment alignment: %s", args.alignment_l)
            AlignmentCleaner.clean_alignment(
                args.alignment_l,
                get_cleaned_output_path(args.alignment_l),
                sequences_to_remove
            )

        # Clean tree files
        if args.tree_s:
            logger.debug("Cleaning S segment tree: %s", args.tree_s)
            TreeCleaner.clean_tree(
                args.tree_s,
                get_cleaned_output_path(args.tree_s),
                sequences_to_remove,
                get_cleaned_output_path(args.alignment_s) if args.alignment_s else None
            )

        if args.tree_l:
            logger.debug("Cleaning L segment tree: %s", args.tree_l)
            TreeCleaner.clean_tree(
                args.tree_l,
                get_cleaned_output_path(args.tree_l),
                sequences_to_remove,
                get_cleaned_output_path(args.alignment_l) if args.alignment_l else None
            )
    else:
        logger.info("No sequences marked for removal from starting alignment")

    logger.info("Processing complete!")
    logger.info("="*50)
    logger.info("SUMMARY:")
    if processor.df is not None:
        logger.info("- Processed %d total sequences from Excel file", len(processor.df))
    logger.info("- Filtered to %d sequences for output", len(filtered_data))
    logger.info("- Generated YAML config: %s", f"{args.output}_filter.yml")
    logger.info("- Generated FASTA files: %s and %s", f"{args.output}_L.fasta", f"{args.output}_S.fasta")
    if sequences_to_remove:
        logger.info("- Removed %d sequences from base alignments/trees", len(sequences_to_remove))
    logger.info("="*50)


if __name__ == "__main__":
    main()