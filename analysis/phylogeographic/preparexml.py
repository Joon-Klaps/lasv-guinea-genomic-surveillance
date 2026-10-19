#!/usr/bin/env python3
"""
XML Generator for BEAST Analysis

This script generates a populated XML file from a Jinja2 template for BEAST analysis.
It reads sequence alignment files, extracts metadata from Excel files, and creates
the necessary XML components for phylogeographic analysis.

Usage:
    python preparexml.py --alignment ALIGNMENT_FILE --metadata METADATA_FILE --output OUTPUT_PREFIX --segment {L,S} [--template TEMPLATE_FILE]
    python preparexml.py --alignment ALIGNMENT_FILE --metadata METADATA_FILE --output OUTPUT_PREFIX --segment {L,S} --codon-partitioned

The script supports two main template types:
- Standard template (xml-templates/template_jinja.xml): Basic phylogeographic analysis
- Codon-partitioned template (xml-templates/template_codon_partitioned.xml): Separate GTR models for each codon position

Requirements:
    - pandas
    - jinja2
    - biopython
    - openpyxl (for Excel reading)
"""

import argparse
import sys
import os
from datetime import datetime
from typing import Optional, Tuple, List, Dict, Any
import pandas as pd
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader
import re


# Minimal ISO 3166 alpha-3 codes for likely countries in LASV datasets
COUNTRY_TO_ISO3 = {
    'guinea': 'GIN',
    'sierra leone': 'SLE',
    'liberia': 'LBR',
    'nigeria': 'NGA',
    'benin': 'BEN',
    'togo': 'TGO',
    'ghana': 'GHA',
    "cote d'ivoire": 'CIV',
    'ivory coast': 'CIV',
    'mali': 'MLI',
    'burkina faso': 'BFA',
    'niger': 'NER',
    'cameroon': 'CMR',
    'central african republic': 'CAF',
    'gabon': 'GAB',
    'equatorial guinea': 'GNQ',
}


def iso3_from_country(country: Any) -> str:
    """Return ISO-3166 alpha-3 code for a country name, or '' if unknown."""
    if pd.isna(country):
        return ''
    name = str(country).strip().lower()
    # normalize common punctuation variants
    name = name.replace('côte d\'ivoire', "cote d'ivoire")
    return COUNTRY_TO_ISO3.get(name, '')


class ExitCodes:
    """Exit codes for different error conditions"""
    SUCCESS = 0
    FILE_NOT_FOUND = 1
    INVALID_FORMAT = 2
    MISSING_DATA = 3
    TEMPLATE_ERROR = 4
    GENERAL_ERROR = 5


class TaxonData:
    """Data class for storing taxon information"""
    def __init__(self, taxon_id: str, date_value: float, latitude: float, longitude: float,
                 uncertainty: Optional[float] = None, needs_tip_sampling: bool = False) -> None:
        self.id = taxon_id
        self.date_value = date_value
        self.latitude = latitude
        self.longitude = longitude
        self.uncertainty = uncertainty
        self.needs_tip_sampling = needs_tip_sampling


class SequenceData:
    """Data class for storing sequence information"""
    def __init__(self, taxon_id: str, sequence: str) -> None:
        self.taxon_id = taxon_id
        self.sequence = sequence


def parse_date_string(date_str: Any) -> Tuple[Optional[float], Optional[float], bool]:
    """
    Parse various date formats and return decimal year and uncertainty.

    Args:
        date_str: Date string in various formats (YYYY, YYYY-MM, YYYY-MM-DD)

    Returns:
        tuple: (decimal_year, uncertainty, needs_tip_sampling)
    """
    if pd.isna(date_str) or str(date_str).strip() == "":
        return None, None, False

    date_str = str(date_str).strip()

    # Pattern for YYYY-MM-DD
    if re.match(r'^\d{4}-\d{2}-\d{2}$', date_str):
        try:
            dt = datetime.strptime(date_str, '%Y-%m-%d')
            day_of_year = dt.timetuple().tm_yday
            decimal_year = dt.year + (day_of_year - 1) / 365.25
            return decimal_year, None, False
        except ValueError:
            return None, None, False

    # Pattern for YYYY-MM
    elif re.match(r'^\d{4}-\d{2}$', date_str):
        try:
            year, month = map(int, date_str.split('-'))
            # Use middle of the month
            dt = datetime(year, month, 15)
            day_of_year = dt.timetuple().tm_yday
            decimal_year = year + (day_of_year - 1) / 365.25
            uncertainty = 1/12  # Monthly uncertainty
            return decimal_year, uncertainty, True
        except ValueError:
            return None, None, True

    # Pattern for YYYY
    elif re.match(r'^\d{4}$', date_str):
        try:
            year = int(date_str)
            decimal_year = year
            uncertainty = 1 # Yearly uncertainty, going forward
            return decimal_year, uncertainty, True
        except ValueError:
            return None, None, True

    return None, None, False


def read_alignment_file(alignment_path: str) -> Tuple[Optional[List[str]], Optional[Dict[str, str]], Optional[int]]:
    """
    Read alignment file and extract sequence IDs and sequences.

    Args:
        alignment_path: Path to FASTA alignment file

    Returns:
        tuple: (sequence_ids, sequences_dict, sequence_length)
    """
    sequences = {}
    sequence_ids = []
    sequence_length = None

    try:
        for record in SeqIO.parse(alignment_path, "fasta"):
            sequences[record.id] = str(record.seq)
            sequence_ids.append(record.id)
            if sequence_length is None:
                sequence_length = len(record.seq)
    except (IOError, ValueError) as e:
        print(f"Error reading alignment file: {e}", file=sys.stderr)
        return None, None, None

    return sequence_ids, sequences, sequence_length


def read_metadata_file(metadata_path: str) -> Optional[pd.DataFrame]:
    """
    Read metadata Excel file and return DataFrame.

    Args:
        metadata_path: Path to Excel metadata file

    Returns:
        pandas.DataFrame or None if error
    """
    try:
        # Try reading the Excel file
        df = pd.read_excel(metadata_path)
        return df
    except (IOError, ValueError, pd.errors.EmptyDataError) as e:
        print(f"Error reading metadata file: {e}", file=sys.stderr)
        return None


def match_sequences_to_metadata(sequence_ids: List[str], metadata_df: pd.DataFrame) -> Optional[Dict[str, Dict[str, Any]]]:
    """
    Match sequence IDs to metadata records.

    Args:
        sequence_ids: List of sequence IDs from alignment
        metadata_df: DataFrame with metadata

    Returns:
        dict: Mapping of sequence_id to metadata row
    """
    matches = {}

    # Common column names that might contain sequence IDs
    id_col = 'alignment_name'
    id_col2 = 'old_alignment_name'
    date_col = 'Collection_Date'
    lat_col = 'latitude'
    lon_col = 'longitude'
    extra_cols = ['Accession', 'Country', 'admin1', 'location', 'Segment', "Host"]

    seq_set = set(sequence_ids)
    meta_set1 = set(metadata_df[id_col].astype(str).str.strip())
    meta_set2 = set(metadata_df[id_col2].astype(str).str.strip())
    meta_set_combined = meta_set1.union(meta_set2)

    missing_ids = seq_set - meta_set_combined

    if missing_ids:
        print(f"WARNING: The following sequence IDs ({len(missing_ids)}) are missing in metadata", file=sys.stderr)
        with open('missing_ids.txt', 'w', encoding='utf-8') as f:
            f.write(f"Missing sequence IDs ({len(missing_ids)}):\n")
            f.write('\n'.join(sorted(missing_ids)))

    # Find matches in both columns
    matching_ids1 = metadata_df[metadata_df[id_col].isin(seq_set)]
    matching_ids2 = metadata_df[metadata_df[id_col2].isin(seq_set)]

    if not matching_ids1.empty:
        df_ids1 = matching_ids1.set_index(id_col)
        df_ids1 = df_ids1.groupby(df_ids1.index).first()
        print(f"Found {len(df_ids1)} matching in {id_col} for metadata records for {len(sequence_ids)} sequence IDs.")
        keep_cols = [c for c in [date_col, lat_col, lon_col] + extra_cols if c in df_ids1.columns]
        matches.update(df_ids1[keep_cols].to_dict('index'))

    if not matching_ids2.empty:
        df_ids2 = matching_ids2.set_index(id_col2)
        df_ids2 = df_ids2.groupby(df_ids2.index).first()
        print(f"Found {len(df_ids2)} matching in {id_col2} for metadata records for {len(sequence_ids)} sequence IDs.")
        keep_cols = [c for c in [date_col, lat_col, lon_col] + extra_cols if c in df_ids2.columns]
        matches.update(df_ids2[keep_cols].to_dict('index'))

    return matches


def generate_xml(template_path: str, taxa_data: List[TaxonData], sequences_data: List[SequenceData],
                 sequence_length: int, output_prefix: str) -> Optional[str]:
    """
    Generate XML file using Jinja2 template.

    Args:
        template_path: Path to Jinja2 template file
        taxa_data: List of TaxonData objects
        sequences_data: List of SequenceData objects
        sequence_length: Length of sequences
        output_prefix: Prefix for output files

    Returns:
        str: Generated XML content or None if error
    """
    try:
        # Set up Jinja2 environment
        template_dir = os.path.dirname(template_path)
        template_name = os.path.basename(template_path)

        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template(template_name)

        # Make sure output_prefix is only basename, without directories
        output_prefix = os.path.basename(output_prefix)

        # Render template
        xml_content = template.render(
            taxa=taxa_data,
            sequences=sequences_data,
            sequence_length=sequence_length,
            output_prefix=output_prefix
        )

        return xml_content

    except (IOError, ValueError) as e:
        print(f"Error generating XML: {e}", file=sys.stderr)
        return None


def main() -> int:
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Generate BEAST XML file from alignment and metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--alignment', '-a',
        required=True,
        help='Path to FASTA alignment file'
    )

    parser.add_argument(
        '--aln-only',
        action='store_true',
        help='Only generate a subset of the algniment'
    )

    parser.add_argument(
        '--metadata', '-m',
        default="../../data/curated-metadata/LASV_GENBANK_20250725_curated.raw.xlsx",
        help='Path to Excel metadata file'
    )

    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output prefix for generated files'
    )

    parser.add_argument(
        '--template', '-t',
        default='xml-templates/template_codon_partitioned.xml',
        help='Path to Jinja2 XML template file'
    )

    parser.add_argument(
        '--segment', '-s',
        choices=['L', 'S'],
        required=True,
        help='Segment code (L or S) used for column naming in exported metadata'
    )

    parser.add_argument(
        '--skip-countries',
        nargs='+',
        default=[],
        help='List of country names to skip sequences from (e.g., NIG CIV LBR)'
    )

    parser.add_argument(
        '--codon-partitioned',
        action='store_true',
        help='Use codon-partitioned template (overrides --template)'
    )

    args = parser.parse_args()

    # Set template based on options
    if args.codon_partitioned:
        args.template = 'xml-templates/template_codon_partitioned.xml'

    # Check if files exist
    if not os.path.exists(args.alignment):
        print(f"Error: Alignment file not found: {args.alignment}", file=sys.stderr)
        return ExitCodes.FILE_NOT_FOUND

    if not os.path.exists(args.metadata):
        print(f"Error: Metadata file not found: {args.metadata}", file=sys.stderr)
        return ExitCodes.FILE_NOT_FOUND

    if not os.path.exists(args.template):
        print(f"Error: Template file not found: {args.template}", file=sys.stderr)
        return ExitCodes.FILE_NOT_FOUND

    # Read alignment file
    print("Reading alignment file...")
    sequence_ids, sequences_dict, sequence_length = read_alignment_file(args.alignment)
    if sequence_ids is None or sequences_dict is None or sequence_length is None:
        return ExitCodes.INVALID_FORMAT

    print(f"Found {len(sequence_ids)} sequences of length {sequence_length}")

    # Read metadata file
    print("Reading metadata file...")
    metadata_df = read_metadata_file(args.metadata)
    if metadata_df is None:
        return ExitCodes.INVALID_FORMAT

    print(f"Found {len(metadata_df)} metadata records")

    # Match sequences to metadata
    print("Matching sequences to metadata...")
    matches = match_sequences_to_metadata(sequence_ids, metadata_df)
    if matches is None:
        return ExitCodes.MISSING_DATA

    print(f"Successfully matched {len(matches)} sequences to metadata")

    if len(matches) == 0:
        print("Warning: No sequences could be matched to metadata", file=sys.stderr)
        return ExitCodes.MISSING_DATA

    # Create taxa data
    taxa_data = []
    sequences_data = []
    skipped_countries = []
    export_rows: List[Dict[str, Any]] = []

    for seq_id in sequence_ids:
        if seq_id in matches:
            metadata = matches[seq_id]

            # Parse date
            date_value, uncertainty, needs_tip_sampling = parse_date_string(metadata['Collection_Date'])

            if date_value is None:
                print(f"Warning: Could not parse date for {seq_id}: {metadata['Collection_Date']}", file=sys.stderr)
                continue

            # Validate coordinates
            try:
                lat = float(metadata['latitude'])
                lon = float(metadata['longitude'])
            except (ValueError, TypeError):
                print(f"Warning: Invalid coordinates for {seq_id}: lat={metadata['latitude']}, lon={metadata['longitude']}",
                      file=sys.stderr)
                continue

            # check if has a country in the first place
            if not metadata.get('Country'):
                print(f"Warning: Missing country for {seq_id}, skipping.", file=sys.stderr)
                continue

            # check if it has latitude and longitude data
            if pd.isna(metadata.get('latitude')) or pd.isna(metadata.get('longitude')):
                print(f"Warning: Missing latitude/longitude for {seq_id}, skipping.", file=sys.stderr)
                continue

            # Skip countries if specified
            if any(country_code.upper() == metadata.get('Country', '').upper() for country_code in args.skip_countries):
                skipped_countries.append(seq_id)
                continue

            # Create taxon data
            taxon = TaxonData(
                taxon_id=seq_id,
                date_value=date_value,
                latitude=lat,
                longitude=lon,
                uncertainty=uncertainty,
                needs_tip_sampling=needs_tip_sampling
            )
            taxa_data.append(taxon)

            # Create sequence data
            sequence = SequenceData(
                taxon_id=seq_id,
                sequence=sequences_dict[seq_id]
            )
            sequences_data.append(sequence)

            # Collect export metadata row (only for included taxa)
            country_val = metadata.get('Country', None)
            admin1_val = metadata.get('admin1', None)
            location_val = metadata.get('location', None)
            accession_val = metadata.get('Accession', None)
            host = metadata.get('Host', None)
            co_val = iso3_from_country(country_val)

            seg_col = f"segment_{args.segment}"
            acc_col = f"accession_{args.segment}"
            export_rows.append({
                # According to prior convention, this column contains the FASTA taxon label
                seg_col: seq_id,
                acc_col: accession_val,
                'Host': host,
                'Collection_Date': metadata['Collection_Date'],
                'CO': co_val,
                'country': country_val,
                'admin1': admin1_val,
                'location': location_val,
                'latitude': lat,
                'longitude': lon,
            })

    print(f"Prepared data for {len(taxa_data)} taxa")
    if skipped_countries:
        print(f"Skipped {len(skipped_countries)} sequences from specified countries due to --skip-countries {args.skip_countries}")

    # Export simple metadata CSV for included taxa
    try:
        if export_rows:
            df_export = pd.DataFrame(export_rows)
            # Ensure all required columns exist (even if empty)
            required_cols = [f'segment_{args.segment}', f'accession_{args.segment}','Collection_Date', 'Host', 'CO', 'country', 'admin1', 'location', 'latitude', 'longitude']
            for col in required_cols:
                if col not in df_export.columns:
                    df_export[col] = ''
            df_export = df_export[required_cols]
            meta_out = f"{args.output}.metadata_{args.segment}.csv"
            df_export.to_csv(meta_out, index=False)
            print(f"Exported metadata for {len(df_export)} taxa -> {meta_out}")
        else:
            print("No taxa available for metadata export.")
    except (OSError, IOError, ValueError) as e:
        print(f"Warning: Failed to export metadata CSV: {e}", file=sys.stderr)

    # Only generate alignment subset if requested
    if args.aln_only:
        aln_out = f"{args.output}.subset_{args.segment}.fasta"
        try:
            with open(aln_out, 'w', encoding='utf-8') as f:
                for seq in sequences_data:
                    f.write(f">{seq.taxon_id}\n{seq.sequence}\n")
            print(f"Exported alignment subset for {len(sequences_data)} taxa -> {aln_out}")
        except (OSError, IOError) as e:
            print(f"Warning: Failed to export alignment subset: {e}", file=sys.stderr)
        return ExitCodes.SUCCESS

    # Generate XML
    print("Generating XML...")
    xml_content = generate_xml(args.template, taxa_data, sequences_data, sequence_length, args.output)
    if xml_content is None:
        return ExitCodes.TEMPLATE_ERROR

    # Write output file
    output_file = f"{args.output}.xml"
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(xml_content)
        print(f"Successfully generated XML file: {output_file}")
    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        return ExitCodes.GENERAL_ERROR

    return ExitCodes.SUCCESS


if __name__ == "__main__":
    sys.exit(main())