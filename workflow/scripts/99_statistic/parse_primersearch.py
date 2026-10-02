import argparse
import re
from typing import Dict, List, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse EMBOSS primersearch output(s) into a contig -> primer_name/primer_match/"
        "amplimer_length/mismatch table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--strict", required=True, help="primersearch output with -mismatchpercent 0")
    parser.add_argument("--relaxed", required=True, help="primersearch output with -mismatchpercent 10")
    parser.add_argument(
        "--output", required=True,
        help="Output TSV: Contig, primer_name, primer_match, amplimer_length, mismatch",
    )
    return parser.parse_args()


def parse_primersearch(path: str) -> List[Tuple[str, str, int, int]]:
    """Extract (contig, primer_name, amplimer_length, mismatch) tuples from a primersearch output file.

    mismatch is the sum of the forward-strand and reverse-strand mismatch counts for that amplimer.
    """
    hits = []
    primer_name = None
    contig = None
    mismatch_sum = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()

            primer_header = re.match(r"Primer name (\S+)", line)
            if primer_header:
                primer_name = primer_header.group(1)
                continue

            seq_match = re.match(r"Sequence:\s*(\S+)", line)
            if seq_match:
                contig = seq_match.group(1)
                mismatch_sum = 0
                continue

            mismatch_match = re.search(r"with (\d+) mismatches?", line)
            if mismatch_match:
                mismatch_sum += int(mismatch_match.group(1))
                continue

            length_match = re.match(r"Amplimer length:\s*(\d+)\s*bp", line)
            if length_match and contig and primer_name:
                hits.append((contig, primer_name, int(length_match.group(1)), mismatch_sum))
                contig = None

    return hits


def build_tag_table(strict_path: str, relaxed_path: str) -> pd.DataFrame:
    """Combine strict/relaxed hits into one row per contig, preferring strict when both match."""
    strict_hits = parse_primersearch(strict_path)
    relaxed_hits = parse_primersearch(relaxed_path)

    contig_hits: Dict[str, Dict[str, List[Tuple[str, int, int]]]] = {}
    for contig, primer, length, mismatch in strict_hits:
        contig_hits.setdefault(contig, {"strict": [], "relaxed": []})["strict"].append((primer, length, mismatch))
    for contig, primer, length, mismatch in relaxed_hits:
        contig_hits.setdefault(contig, {"strict": [], "relaxed": []})["relaxed"].append((primer, length, mismatch))

    rows = []
    for contig, levels in contig_hits.items():
        primer_match = "strict" if levels["strict"] else "relaxed"
        chosen = levels["strict"] if levels["strict"] else levels["relaxed"]
        chosen = sorted(chosen, key=lambda x: x[0])

        primer_name = ";".join(p for p, _, _ in chosen)
        amplimer_length = ";".join(str(l) for _, l, _ in chosen)
        mismatch = ";".join(str(m) for _, _, m in chosen)
        rows.append((contig, primer_name, primer_match, amplimer_length, mismatch))

    return pd.DataFrame(rows, columns=["Contig", "primer_name", "primer_match", "amplimer_length", "mismatch"])


def main() -> None:
    """Parse primersearch outputs and write the contig primer tag table."""
    args = parse_args()
    table = build_tag_table(args.strict, args.relaxed)
    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
