import argparse
from typing import Dict, List, Tuple

import pandas as pd

COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse BLASTN single-primer outfmt 6 results into a contig -> "
        "single_primer_name/single_primer_match/single_primer_mismatch table. A hit is "
        "kept only if it covers the full primer length (no gaps, length == qlen).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--blast", required=True, help="BLASTN outfmt 6 result file (see COLUMNS)")
    parser.add_argument(
        "--output", required=True,
        help="Output TSV: Contig, single_primer_name, single_primer_match, single_primer_mismatch",
    )
    return parser.parse_args()


def load_hits(path: str) -> pd.DataFrame:
    """Load full-length, gapless BLASTN hits (one row per primer x contig alignment)."""
    df = pd.read_csv(path, sep="\t", names=COLUMNS)
    df = df[(df["gapopen"] == 0) & (df["length"] == df["qlen"])]
    return df


def build_tag_table(blast_path: str) -> pd.DataFrame:
    """Combine hits into one row per contig; strict = 0 mismatches, relaxed = any mismatches."""
    hits = load_hits(blast_path)

    contig_hits: Dict[str, Dict[str, List[Tuple[str, int]]]] = {}
    for _, row in hits.iterrows():
        contig = row["sseqid"]
        primer = row["qseqid"]
        mismatch = int(row["mismatch"])
        level = "strict" if mismatch == 0 else "relaxed"
        contig_hits.setdefault(contig, {"strict": [], "relaxed": []})[level].append((primer, mismatch))

    rows = []
    for contig, levels in contig_hits.items():
        primer_match = "strict" if levels["strict"] else "relaxed"
        chosen = levels["strict"] if levels["strict"] else levels["relaxed"]
        chosen = sorted(set(chosen), key=lambda x: x[0])

        primer_name = ";".join(p for p, _ in chosen)
        mismatch = ";".join(str(m) for _, m in chosen)
        rows.append((contig, primer_name, primer_match, mismatch))

    return pd.DataFrame(rows, columns=["Contig", "single_primer_name", "single_primer_match", "single_primer_mismatch"])


def main() -> None:
    """Parse BLASTN single-primer results and write the contig single-primer tag table."""
    args = parse_args()
    table = build_tag_table(args.blast)
    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
