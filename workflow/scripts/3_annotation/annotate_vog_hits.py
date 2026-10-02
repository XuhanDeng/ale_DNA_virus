import argparse
from typing import Dict

import pandas as pd

DIAMOND_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate DIAMOND blastp hits against the VOG protein database with VOG "
        "group membership, functional annotation, and LCA taxonomy; keep only the best hit "
        "(highest bitscore, ties broken by lowest e-value) per query protein",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--diamond", required=True, help="DIAMOND blastp outfmt 6 results file")
    parser.add_argument("--members", required=True, help="vog.members.tsv")
    parser.add_argument("--annotations", required=True, help="vog.annotations.tsv")
    parser.add_argument("--lca", required=True, help="vog.lca.tsv")
    parser.add_argument("--output", required=True, help="Annotated output TSV")
    return parser.parse_args()


def load_protein_to_group(members_path: str) -> Dict[str, str]:
    """Build a protein ID -> VOG group name lookup from vog.members.tsv."""
    members = pd.read_csv(members_path, sep="\t")
    protein_to_group: Dict[str, str] = {}
    for group, protein_ids in zip(members["#GroupName"], members["ProteinIDs"]):
        for protein_id in protein_ids.split(","):
            protein_to_group[protein_id] = group
    return protein_to_group


def load_best_hits(diamond_path: str) -> pd.DataFrame:
    """Load DIAMOND hits and keep only the best hit per query (highest bitscore, then lowest e-value)."""
    hits = pd.read_csv(diamond_path, sep="\t", names=DIAMOND_COLUMNS)
    hits = hits.sort_values(["bitscore", "evalue"], ascending=[False, True])
    return hits.drop_duplicates("qseqid", keep="first")


def main() -> None:
    """Annotate best-hit DIAMOND results with VOG group, function, and LCA taxonomy."""
    args = parse_args()

    best_hits = load_best_hits(args.diamond)
    protein_to_group = load_protein_to_group(args.members)
    best_hits["vog_group"] = best_hits["sseqid"].map(protein_to_group)

    annotations = pd.read_csv(args.annotations, sep="\t")
    annotated = pd.merge(best_hits, annotations, how="left", left_on="vog_group", right_on="#GroupName")

    lca = pd.read_csv(args.lca, sep="\t")
    annotated = pd.merge(annotated, lca, how="left", left_on="vog_group", right_on="#GroupName", suffixes=("", "_lca"))

    annotated = annotated.drop(columns=[c for c in annotated.columns if c.startswith("#GroupName")])
    annotated.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
