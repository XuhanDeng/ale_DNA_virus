import argparse

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Select VOG group names whose ConsensusFunctionalDescription contains "
        "the keyword phrase 'terminase large subunit' (case-insensitive substring match)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--annotations", required=True, help="vog.annotations.tsv")
    parser.add_argument("--keyword", default="terminase large subunit", help="Keyword phrase to match")
    parser.add_argument("--output", required=True, help="Output list of matching VOG group names")
    return parser.parse_args()


def select_matching_groups(annotations_path: str, keyword: str) -> pd.Series:
    """Return VOG group names whose functional description contains the keyword phrase."""
    annotations = pd.read_csv(annotations_path, sep="\t")
    mask = annotations["ConsensusFunctionalDescription"].str.contains(keyword, case=False, na=False)
    return annotations.loc[mask, "#GroupName"]


def main() -> None:
    """Select and write VOG group names matching the keyword phrase."""
    args = parse_args()
    matching_groups = select_matching_groups(args.annotations, args.keyword)
    matching_groups.to_csv(args.output, index=False, header=False)


if __name__ == "__main__":
    main()
