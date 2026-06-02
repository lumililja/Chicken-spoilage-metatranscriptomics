"""Build ko_all_reactions_4.csv or ko_all_reactions_6.csv from HUMAnN outputs.

This script does the following:
1. Regroup gene families to KO identifiers.
2. Run HUMAnN with the KO pathway/module database.
3. Clean, merge, and annotate the stratified reaction abundance tables.

The external ``humann_regroup_table`` and ``humann`` commands must be installed
and available on PATH for the first two stages.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class TemperatureConfig:
    sample_ids: tuple[str, ...]
    gene_families_template: str
    ko_dir: str
    humann_rebuild_dir: str
    cleaned_reactions_dir: str
    output_csv: str


CONFIGS = {
    4: TemperatureConfig(
        sample_ids=tuple(f"A{i:02d}" for i in range(1, 31)),
        gene_families_template=(
            "humann_output_4_degrees/"
            "{sample_id}_genefamilies.tsv"
        ),
        ko_dir="KO_4",
        humann_rebuild_dir="humann_rebuild_4",
        cleaned_reactions_dir="KO_reactions_pathabundances_4",
        output_csv="ko_all_reactions_4.csv",
    ),
    6: TemperatureConfig(
        sample_ids=tuple(f"A{i:02d}" for i in range(5, 20)),
        gene_families_template=(
            "humann_output_6_degrees/"
            "concat_{sample_id}_genefamilies.tsv"
        ),
        ko_dir="KO_6",
        humann_rebuild_dir="humann_rebuild_6",
        cleaned_reactions_dir="KO_reactions_pathabundances_6",
        output_csv="ko_all_reactions_6.csv",
    ),
}

def require_executable(name: str) -> None:
    """Raise a useful error before starting a stage with a missing tool."""
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required executable '{name}' was not found on PATH. "
            "Activate the environment containing HUMAnN and run the command again."
        )


def run_command(command: list[str], expected_output: Path) -> None:
    """Run a checked external command and confirm that it created its output."""
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)
    if not expected_output.exists():
        raise RuntimeError(f"Command completed but did not create {expected_output}")


def regroup_gene_families(workspace: Path, config: TemperatureConfig) -> None:
    """Create one KO regrouped gene-family table per sample."""
    require_executable("humann_regroup_table")
    ko_dir = workspace / config.ko_dir
    ko_dir.mkdir(parents=True, exist_ok=True)
    custom_mapping = workspace / "utility_mapping/map_ko-expanded_uniref90.txt.gz"

    for sample_id in config.sample_ids:
        input_tsv = workspace / config.gene_families_template.format(
            sample_id=sample_id
        )
        output_tsv = ko_dir / f"ko_{sample_id}.tsv"
        if not input_tsv.exists():
            raise FileNotFoundError(f"Missing gene-family input: {input_tsv}")
        run_command(
            [
                "humann_regroup_table",
                "--input",
                str(input_tsv),
                "--ungrouped",
                "N",
                "--protected",
                "N",
                "--custom",
                str(custom_mapping),
                "--output",
                str(output_tsv),
            ],
            output_tsv,
        )


def rebuild_humann_pathways(workspace: Path, config: TemperatureConfig) -> None:
    """Run HUMAnN on regrouped KO tables to generate reaction abundances."""
    require_executable("humann")
    ko_dir = workspace / config.ko_dir
    rebuild_dir = workspace / config.humann_rebuild_dir
    rebuild_dir.mkdir(parents=True, exist_ok=True)
    modules = workspace / "utility_mapping/ko_pathway_modules.txt"

    for sample_id in config.sample_ids:
        input_tsv = ko_dir / f"ko_{sample_id}.tsv"
        output_tsv = rebuild_dir / f"ko_{sample_id}_3_reactions.tsv"
        if not input_tsv.exists():
            raise FileNotFoundError(f"Missing KO input: {input_tsv}")
        run_command(
            [
                "humann",
                "--input",
                str(input_tsv),
                "--input-format",
                "genetable",
                "--log-level",
                "DEBUG",
                "--pathways-database",
                f"{modules},{modules}",
                "--output",
                str(rebuild_dir),
                "--minpath",
                "off",
                "--gap-fill",
                "off",
            ],
            output_tsv,
        )


def clean_reaction_table(input_tsv: Path, sample_id: str) -> pd.DataFrame:
    """Keep stratified biological reactions and standardize the sample column."""
    df = pd.read_csv(input_tsv, sep="\t", index_col=0)
    if len(df.columns) != 1:
        raise ValueError(f"Expected one abundance column in {input_tsv}, found {df.columns}")
    df = df.loc[~df.index.str.contains("UNMAPPED|UNINTEGRATED|UNGROUPED")]
    df = df.loc[df.index.str.contains(r"\|")]
    return df.rename(columns={df.columns[0]: sample_id})


def read_mapping(path: Path) -> dict[str, str]:
    """Read a colon-delimited mapping while allowing colons inside values."""
    mapping = {}
    with path.open() as mapping_file:
        next(mapping_file, None)
        for line in mapping_file:
            key, value = line.rstrip("\n").split(":", 1)
            mapping[key] = value
    return mapping


def annotate_reactions(workspace: Path, reactions: pd.DataFrame) -> pd.DataFrame:
    """Add reaction ID, bacterium, reaction name, and two-level pathway class."""
    names = read_mapping(workspace / "utility_mapping/all_pathways_names.txt")
    classes = read_mapping(workspace / "utility_mapping/all_pathways_class.txt")

    annotations = reactions.index.to_series().str.split("|", n=1, expand=True)
    reactions["reaction"] = annotations[0].values
    reactions["bacterium"] = annotations[1].values
    reactions["name"] = reactions["reaction"].map(names)

    def format_class(reaction: str) -> str:
        levels = classes[reaction].split(";")
        if len(levels) < 3:
            raise ValueError(f"Expected at least three class levels for {reaction}")
        return f"{levels[1].strip()}, {levels[2].strip()}"

    missing_names = sorted(set(reactions["reaction"]) - names.keys())
    missing_classes = sorted(set(reactions["reaction"]) - classes.keys())
    if missing_names or missing_classes:
        raise KeyError(
            "Missing reaction metadata: "
            f"names={missing_names[:5]}, classes={missing_classes[:5]}"
        )
    reactions["class"] = reactions["reaction"].map(format_class)
    return reactions


def merge_reactions(workspace: Path, config: TemperatureConfig) -> Path:
    """Clean and merge per-sample reaction TSV files into the final CSV."""
    rebuild_dir = workspace / config.humann_rebuild_dir
    cleaned_dir = workspace / config.cleaned_reactions_dir
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    tables = []

    for sample_id in config.sample_ids:
        input_tsv = rebuild_dir / f"ko_{sample_id}_3_reactions.tsv"
        if not input_tsv.exists():
            raise FileNotFoundError(f"Missing HUMAnN reaction table: {input_tsv}")
        table = clean_reaction_table(input_tsv, sample_id)
        table.to_csv(cleaned_dir / f"{sample_id}_reactions.tsv", sep="\t")
        tables.append(table)

    reactions = pd.concat(tables, axis=1).fillna(0)
    reactions = annotate_reactions(workspace, reactions)
    output_csv = workspace / config.output_csv
    reactions.to_csv(output_csv)
    print(f"Wrote {output_csv}")
    return output_csv


def build_ko_all_reactions(
    temperature: int,
    *,
    workspace: str | Path = ".",
) -> Path:
    """Run the configured pipeline and return the final CSV path."""
    try:
        config = CONFIGS[temperature]
    except KeyError as exc:
        raise ValueError(f"Unsupported temperature: {temperature}. Use 4 or 6.") from exc
    workspace = Path(workspace).resolve()
    regroup_gene_families(workspace, config)
    rebuild_humann_pathways(workspace, config)
    return merge_reactions(workspace, config)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--temperature",
        type=int,
        choices=sorted(CONFIGS),
        required=True,
        help="Storage temperature in degrees Celsius.",
    )
    parser.add_argument(
        "--workspace",
        default=".",
        help="Directory containing the HUMAnN inputs and utility_mapping directory.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_ko_all_reactions(
        args.temperature,
        workspace=args.workspace,
    )
