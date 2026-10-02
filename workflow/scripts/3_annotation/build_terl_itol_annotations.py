import argparse
import colorsys
from typing import Dict, List, Set

import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build iTOL DATASET_COLORSTRIP annotation files for the TerL tree's VOG "
        "reference leaves: one ring colored by VOG group, one by LCA taxonomy. Using "
        "colorstrip rings (rather than TREE_COLORS) lets both layers be shown at once "
        "without overlapping each other or the tree's branch colors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--references", required=True, help="terl_vog_references.faa (pre-filter VOG reference FASTA)")
    parser.add_argument("--alignment-input", required=True, help="all_TerL_filtered.faa (sequences that entered the tree)")
    parser.add_argument("--members", required=True, help="vog.members.tsv")
    parser.add_argument("--lca", required=True, help="vog.lca.tsv")
    parser.add_argument("--vog-output", required=True, help="Output iTOL DATASET_COLORSTRIP file for VOG group")
    parser.add_argument("--lca-output", required=True, help="Output iTOL DATASET_COLORSTRIP file for LCA taxonomy")
    return parser.parse_args()


def load_tree_reference_ids(references_path: str, alignment_input_path: str) -> Set[str]:
    """IDs that are both VOG references and survived filtering into the tree."""
    reference_ids = {record.id for record in SeqIO.parse(references_path, "fasta")}
    tree_ids = {record.id for record in SeqIO.parse(alignment_input_path, "fasta")}
    return reference_ids & tree_ids


def load_protein_to_group(members_path: str) -> Dict[str, str]:
    """Build a protein ID -> VOG group name lookup from vog.members.tsv."""
    members = pd.read_csv(members_path, sep="\t")
    protein_to_group: Dict[str, str] = {}
    for group, protein_ids in zip(members["#GroupName"], members["ProteinIDs"]):
        for protein_id in protein_ids.split(","):
            protein_to_group[protein_id] = group
    return protein_to_group


def load_group_to_lca(lca_path: str) -> Dict[str, str]:
    """Build a VOG group name -> LastCommonAncestor_Name lookup from vog.lca.tsv."""
    lca = pd.read_csv(lca_path, sep="\t")
    return dict(zip(lca["#GroupName"], lca["LastCommonAncestor_Name"]))


def assign_colors(categories: List[str]) -> Dict[str, str]:
    """Assign a distinct hex color to each category using evenly spaced hues."""
    unique = sorted(set(categories))
    n = len(unique)
    colors = {}
    for i, category in enumerate(unique):
        hue = i / n if n else 0
        r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.85)
        colors[category] = "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
    return colors


def write_colorstrip(leaf_to_value: Dict[str, str], label: str, output_path: str) -> None:
    """Write an iTOL DATASET_COLORSTRIP annotation file, one ring block color per distinct value.

    Colorstrip rings render outside the tree (next to leaf labels) rather than coloring
    branches directly, so multiple annotation layers can be shown together without overlap.
    """
    color_map = assign_colors(list(leaf_to_value.values()))
    with open(output_path, "w") as handle:
        handle.write("DATASET_COLORSTRIP\n")
        handle.write("SEPARATOR TAB\n")
        handle.write(f"DATASET_LABEL\t{label}\n")
        handle.write("COLOR\t#000000\n")
        handle.write("COLOR_BRANCHES\t0\n")
        handle.write("LEGEND_TITLE\t" + label + "\n")
        legend_values = sorted(color_map)
        handle.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in legend_values) + "\n")
        handle.write("LEGEND_COLORS\t" + "\t".join(color_map[v] for v in legend_values) + "\n")
        handle.write("LEGEND_LABELS\t" + "\t".join(legend_values) + "\n")
        handle.write("DATA\n")
        for leaf, value in leaf_to_value.items():
            handle.write(f"{leaf}\t{color_map[value]}\t{value}\n")


def main() -> None:
    """Build VOG-group and LCA-taxonomy iTOL colorstrip ring files for the tree's reference leaves."""
    args = parse_args()

    tree_reference_ids = load_tree_reference_ids(args.references, args.alignment_input)
    protein_to_group = load_protein_to_group(args.members)
    group_to_lca = load_group_to_lca(args.lca)

    leaf_to_group: Dict[str, str] = {}
    leaf_to_lca: Dict[str, str] = {}
    for leaf_id in tree_reference_ids:
        group = protein_to_group.get(leaf_id)
        if group is None:
            continue
        leaf_to_group[leaf_id] = group
        leaf_to_lca[leaf_id] = group_to_lca.get(group, "Unknown")

    write_colorstrip(leaf_to_group, "VOG group", args.vog_output)
    write_colorstrip(leaf_to_lca, "LCA taxonomy", args.lca_output)


if __name__ == "__main__":
    main()
