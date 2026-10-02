import argparse
from typing import Dict, Tuple

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert a FASTA of >NAME_F/>NAME_R primer records into EMBOSS "
        "primersearch's -infile format (one 'name forward reverse' triplet per line)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta", required=True, help="Primer FASTA with >{pair}_F / >{pair}_R records")
    parser.add_argument("--output", required=True, help="Output primersearch -infile path")
    return parser.parse_args()


def load_primer_pairs(path: str) -> Dict[str, Dict[str, str]]:
    """Group primer sequences by pair name and direction (F/R)."""
    pairs: Dict[str, Dict[str, str]] = {}
    for record in SeqIO.parse(path, "fasta"):
        name, direction = record.id.rsplit("_", 1)
        pairs.setdefault(name, {})[direction] = str(record.seq)
    return pairs


def write_primersearch_infile(pairs: Dict[str, Dict[str, str]], output: str) -> None:
    """Write pairs with both F and R primers to the primersearch -infile format."""
    with open(output, "w") as handle:
        for name, seqs in sorted(pairs.items()):
            if "F" in seqs and "R" in seqs:
                handle.write(f"{name}\t{seqs['F']}\t{seqs['R']}\n")


def main() -> None:
    """Convert primer FASTA to EMBOSS primersearch -infile format."""
    args = parse_args()
    pairs = load_primer_pairs(args.fasta)
    write_primersearch_infile(pairs, args.output)


if __name__ == "__main__":
    main()
