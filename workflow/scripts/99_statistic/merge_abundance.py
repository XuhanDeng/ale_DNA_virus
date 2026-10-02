import argparse
from typing import List, Set

import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge per-sample CoverM abundance tables (TPM/mean/count) into one combined "
        "table, annotated with contig length, virus membership, and primer-detection tags",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Per-sample abundance TSV files")
    parser.add_argument("--samples", nargs="+", required=True, help="Sample names, matching --inputs order")
    parser.add_argument("--length", required=True, help="Contig length TSV (Contig, length)")
    parser.add_argument("--virus-fasta", required=True, help="Identified-viral cluster representative FASTA")
    parser.add_argument(
        "--primer-tags", required=True,
        help="Primer tag TSV (Contig, primer_name, primer_match, amplimer_length, mismatch)",
    )
    parser.add_argument(
        "--single-primer-tags", required=True,
        help="Single-primer tag TSV (Contig, single_primer_name, single_primer_match, single_primer_mismatch)",
    )
    parser.add_argument("--output-csv", required=True, help="Merged output CSV path")
    parser.add_argument("--output-xlsx", required=True, help="Merged output XLSX path")
    return parser.parse_args()


def merge_abundance(paths: List[str], samples: List[str]) -> pd.DataFrame:
    """Outer-join per-sample CoverM tables on contig ID, filling missing values with 0."""
    merged = None
    for path, sample in zip(paths, samples):
        df = pd.read_csv(path, sep="\t")
        df.columns = ["Contig", sample]
        merged = df if merged is None else pd.merge(merged, df, how="outer", on="Contig")
    merged[samples] = merged[samples].fillna(0)
    return merged


def load_virus_ids(fasta_path: str) -> Set[str]:
    """Collect contig IDs present in the identified-viral cluster FASTA."""
    return {record.id for record in SeqIO.parse(fasta_path, "fasta")}


def main() -> None:
    """Merge, annotate, and write per-sample abundance tables as CSV and XLSX."""
    args = parse_args()
    merged = merge_abundance(args.inputs, args.samples)

    length_df = pd.read_csv(args.length, sep="\t")
    merged = pd.merge(merged, length_df, how="left", on="Contig")

    virus_ids = load_virus_ids(args.virus_fasta)
    merged["virus"] = merged["Contig"].isin(virus_ids).map({True: "yes", False: "no"})

    primer_df = pd.read_csv(args.primer_tags, sep="\t")
    merged = pd.merge(merged, primer_df, how="left", on="Contig")

    single_primer_df = pd.read_csv(args.single_primer_tags, sep="\t")
    merged = pd.merge(merged, single_primer_df, how="left", on="Contig")

    cols = [
        "Contig", "length", "virus",
        "primer_name", "primer_match", "amplimer_length", "mismatch",
        "single_primer_name", "single_primer_match", "single_primer_mismatch",
    ] + args.samples
    merged = merged[cols]

    merged.to_csv(args.output_csv, index=False)
    merged.to_excel(args.output_xlsx, index=False)


if __name__ == "__main__":
    main()
