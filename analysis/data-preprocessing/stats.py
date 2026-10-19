#!/usr/bin/env python3

"""
Pairwise sequence comparison utilities and CLI.

Features
--------
- Read alignment FASTA (or unaligned FASTA) and compare selected samples pairwise.
- Select samples by a names file (one label per line), matching labels in FASTA headers.
- Upper-triangle comparisons only: pairs (i, j) where i < j.
- Optional on-the-fly pairwise alignment via pyFAMSA when --unaligned is passed.
- Importable functions for notebook workflows (extract -> align -> compute stats).

Example CLI
-----------
python stats.py --aln S-G0959.fna --names names.txt --out out.tsv
python stats.py --aln S-unaligned.fasta --names names.txt --unaligned --out out.tsv
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence as TypingSequence, Tuple


DNA_VALID = {"A", "C", "G", "T", "N"}


def _open_text(path: Path):
	return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path, "rt", encoding="utf-8")


def read_fasta(path: Path) -> List[Tuple[str, str]]:
	"""Read all FASTA records as (header, uppercased_seq)."""
	records: List[Tuple[str, str]] = []
	header: str | None = None
	chunks: List[str] = []

	try:
		with _open_text(path) as handle:
			for raw in handle:
				line = raw.strip()
				if not line:
					continue
				if line.startswith(">"):
					if header is not None:
						records.append((header, "".join(chunks).upper()))
					header = line[1:].strip()
					chunks = []
				else:
					chunks.append(line)
			if header is not None:
				records.append((header, "".join(chunks).upper()))
	except FileNotFoundError as exc:
		raise SystemExit(f"FASTA not found: {path}") from exc

	if not records:
		raise SystemExit(f"No sequences found in FASTA: {path}")
	return records


def read_names(path: Path) -> List[str]:
	"""Read one sample label per line from names file."""
	try:
		with path.open("r", encoding="utf-8") as handle:
			names = [line.strip() for line in handle if line.strip()]
	except FileNotFoundError as exc:
		raise SystemExit(f"Names file not found: {path}") from exc

	if not names:
		raise SystemExit(f"No names found in: {path}")
	return names


def filter_records_by_names(records: TypingSequence[Tuple[str, str]], names: TypingSequence[str]) -> List[Tuple[str, str]]:
	"""Return records whose headers contain any requested name.

	Matching is substring-based and keeps input FASTA order.
	"""
	selected: List[Tuple[str, str]] = []
	for header, seq in records:
		if any(name in header for name in names):
			selected.append((header, seq))
	return selected


def _short_name_from_header(header: str, names: TypingSequence[str]) -> str:
	"""Return first matching requested label; fallback to first token of header."""
	for name in names:
		if name in header:
			return name
	return header.split()[0]


def ensure_equal_length(records: TypingSequence[Tuple[str, str]]) -> None:
	"""Validate records are an alignment (equal lengths)."""
	lengths = {len(seq) for _, seq in records}
	if len(lengths) > 1:
		raise SystemExit(
			"Input appears unaligned (sequence lengths differ). "
			"Pass --unaligned to align each pair with pyFAMSA first."
		)


def align_pair_famsa(seq1_id: str, seq1: str, seq2_id: str, seq2: str) -> Tuple[str, str]:
	"""Align two sequences with pyFAMSA and return aligned strings.

	Keeps pyFAMSA import local to avoid requiring it unless alignment is requested.
	"""
	try:
		from pyfamsa import Aligner, Sequence  # type: ignore
	except Exception as exc:
		raise SystemExit(
			"pyFAMSA is required for --unaligned mode. Install with: pip install pyfamsa"
		) from exc

	aligner = Aligner()
	result = aligner.align([
		Sequence(seq1_id.encode("utf-8"), seq1.encode("utf-8")),
		Sequence(seq2_id.encode("utf-8"), seq2.encode("utf-8")),
	])

	aligned: Dict[str, str] = {}
	for item in result:
		rid = item.id.decode("utf-8")
		rseq = item.sequence.decode("utf-8").upper()
		aligned[rid] = rseq

	return aligned[seq1_id], aligned[seq2_id]


def compute_stats(seq1: str, seq2: str) -> Dict[str, int | float | None]:
	"""Compute pairwise comparison stats on aligned sequences."""
	if len(seq1) != len(seq2):
		raise ValueError("Aligned sequences must have equal lengths")

	aln_len = len(seq1)
	comparable = 0
	matches = 0
	mismatches = 0

	for a, b in zip(seq1, seq2):
		if a in DNA_VALID and b in DNA_VALID:
			comparable += 1
			if a == b and a != "N":
				matches += 1
			elif a != b and a != "N" and b != "N":
				mismatches += 1

	mismatch_rate = (mismatches / comparable) if comparable else None

	return {
		"alignment_length": aln_len,
		"comparable_sites": comparable,
		"matches": matches,
		"mismatches": mismatches,
		"mismatch_rate": mismatch_rate,
		"seq1_length": sum(seq1.count(x) for x in "ACGTN"),
		"seq2_length": sum(seq2.count(x) for x in "ACGTN"),
		"seq1_gaps": seq1.count("-"),
		"seq2_gaps": seq2.count("-"),
		"seq1_Ns": seq1.count("N"),
		"seq2_Ns": seq2.count("N"),
	}


def pairwise_stats(
	records: TypingSequence[Tuple[str, str]],
	names: TypingSequence[str],
	*,
	unaligned: bool = False,
) -> List[List[str]]:
	"""Compute upper-triangle pairwise stats rows for selected records."""
	selected = filter_records_by_names(records, names)
	if len(selected) < 2:
		raise SystemExit(
			f"Need at least 2 matched sequences for pairwise comparison; found {len(selected)}"
		)

	if not unaligned:
		ensure_equal_length(selected)

	rows: List[List[str]] = []
	for (h1, s1), (h2, s2) in itertools.combinations(selected, 2):
		id1 = _short_name_from_header(h1, names)
		id2 = _short_name_from_header(h2, names)

		if unaligned:
			a1, a2 = align_pair_famsa(id1, s1.replace("-", ""), id2, s2.replace("-", ""))
			stats = compute_stats(a1, a2)
		else:
			stats = compute_stats(s1, s2)

		rows.append([
			id1,
			id2,
			str(stats["alignment_length"]),
			str(stats["comparable_sites"]),
			str(stats["matches"]),
			str(stats["mismatches"]),
			f"{stats['mismatch_rate']:.6f}" if stats["mismatch_rate"] is not None else "NA",
			str(stats["seq1_length"]),
			str(stats["seq2_length"]),
			str(stats["seq1_gaps"]),
			str(stats["seq2_gaps"]),
			str(stats["seq1_Ns"]),
			str(stats["seq2_Ns"]),
		])
	return rows


def tsv_header() -> List[str]:
	return [
		"seq1_id",
		"seq2_id",
		"alignment_length",
		"comparable_sites",
		"matches",
		"mismatches",
		"mismatch_rate",
		"seq1_length",
		"seq2_length",
		"seq1_gaps",
		"seq2_gaps",
		"seq1_Ns",
		"seq2_Ns",
	]


def write_tsv(out_path: Path | None, header: TypingSequence[str], rows: Iterable[TypingSequence[str]]) -> None:
	line_header = "\t".join(header) + "\n"
	if out_path:
		out_path.parent.mkdir(parents=True, exist_ok=True)
		with out_path.open("w", encoding="utf-8") as handle:
			handle.write(line_header)
			for row in rows:
				handle.write("\t".join(row) + "\n")
	else:
		sys.stdout.write(line_header)
		for row in rows:
			sys.stdout.write("\t".join(row) + "\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	p.add_argument("--aln", required=True, help="Path to FASTA (.fa/.fasta[.gz])")
	p.add_argument("--names", required=True, help="Path to names.txt (one label per line)")
	p.add_argument("--out", default=None, help="Output TSV path; default STDOUT")
	p.add_argument(
		"--unaligned",
		action="store_true",
		help="Input sequences are unaligned; align each pair with pyFAMSA before stats",
	)
	return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
	args = parse_args(argv)
	aln_path = Path(args.aln)
	names_path = Path(args.names)
	out_path = Path(args.out) if args.out else None

	records = read_fasta(aln_path)
	names = read_names(names_path)
	rows = pairwise_stats(records, names, unaligned=args.unaligned)
	write_tsv(out_path, tsv_header(), rows)


if __name__ == "__main__":
	main()
