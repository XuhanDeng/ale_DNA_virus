import argparse

from Bio import SeqIO

HYPOTHETICAL_PRODUCT = "hypothetical protein"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Remove GenBank records that crash pharokka multiplot (IndexError in "
        "np.quantile on an empty length array), which happens on short/sparse contigs even "
        "when they have a non-hypothetical CDS. Keeps only records that are both >= "
        "--min-length and have at least one non-hypothetical CDS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--genbank", required=True, help="Pharokka pharokka.gbk input file")
    parser.add_argument("--output", required=True, help="Filtered GenBank output file")
    parser.add_argument("--min-length", type=int, default=10000, help="Minimum contig length (bp)")
    return parser.parse_args()


def has_annotated_cds(record) -> bool:
    """True if the record has at least one CDS whose product is not 'hypothetical protein'."""
    for feature in record.features:
        if feature.type != "CDS":
            continue
        product = feature.qualifiers.get("product", [""])[0]
        if product != HYPOTHETICAL_PRODUCT:
            return True
    return False


def filter_records(genbank_path: str, output_path: str, min_length: int) -> None:
    """Write only GenBank records that are >= min_length and have a non-hypothetical CDS."""
    kept, dropped = 0, []
    with open(output_path, "w") as handle:
        for record in SeqIO.parse(genbank_path, "genbank"):
            if len(record.seq) >= min_length and has_annotated_cds(record):
                SeqIO.write(record, handle, "genbank")
                kept += 1
            else:
                dropped.append(record.id)

    print(f"Kept {kept} record(s) >= {min_length}bp with >=1 non-hypothetical CDS.")
    if dropped:
        print(f"Dropped {len(dropped)} record(s): {', '.join(dropped)}")


def main() -> None:
    """Filter records that would crash pharokka multiplot out of a GenBank file."""
    args = parse_args()
    filter_records(args.genbank, args.output, args.min_length)


if __name__ == "__main__":
    main()
