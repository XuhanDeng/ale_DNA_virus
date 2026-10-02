import argparse
import os

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Split a GenBank file into N roughly-equal chunk files by record count, "
        "so pharokka multiplot (which has no multi-threading option and processes one "
        "GenBank file's records sequentially) can be run on each chunk in parallel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--genbank", required=True, help="Input GenBank file to split")
    parser.add_argument("--chunks", type=int, required=True, help="Number of chunk files to create")
    parser.add_argument("--outdir", required=True, help="Directory to write chunk_0.gbk .. chunk_{N-1}.gbk")
    return parser.parse_args()


def split_into_chunks(genbank_path: str, chunks: int, outdir: str) -> None:
    """Round-robin assign records to N chunk directories, one record per file
    (record_<n>.gbk) so a single crashing contig only costs that one plot."""
    os.makedirs(outdir, exist_ok=True)
    chunk_dirs = [os.path.join(outdir, f"chunk_{i}") for i in range(chunks)]
    for chunk_dir in chunk_dirs:
        os.makedirs(chunk_dir, exist_ok=True)

    counts = [0] * chunks
    for i, record in enumerate(SeqIO.parse(genbank_path, "genbank")):
        chunk = i % chunks
        record_path = os.path.join(chunk_dirs[chunk], f"record_{counts[chunk]}.gbk")
        with open(record_path, "w") as handle:
            SeqIO.write(record, handle, "genbank")
        counts[chunk] += 1


def main() -> None:
    """Split a GenBank file into N chunk files."""
    args = parse_args()
    split_into_chunks(args.genbank, args.chunks, args.outdir)


if __name__ == "__main__":
    main()
