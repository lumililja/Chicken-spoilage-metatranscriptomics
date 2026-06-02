"""Create the final a)-d) stacked pathway abundance plots for chicken samples."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


COLOR_MAP = {
    "Central carbohydrate metabolism": "#1f77b4",
    "Other carbohydrate metabolism": "#aec7e8",
    "Arginine and proline metabolism": "#ff7f0e",
    "Lysine metabolism": "#ffbb78",
    "Cysteine and methionine metabolism": "#2ca02c",
    "Serine and threonine metabolism": "#98df8a",
    "Branched-chain amino acid metabolism": "#d62728",
    "Aromatic amino acid metabolism": "#ff9896",
    "Histidine metabolism": "#9467bd",
    "Other amino acid metabolism": "#c5b0d5",
    "Biosynthesis of phytochemical compounds": "#8c564b",
    "Purine metabolism": "#c49c94",
    "Pyrimidine metabolism": "#e377c2",
    "Lipopolysaccharide metabolism": "#f7b6d2",
    "Fatty acid metabolism": "#7f7f7f",
    "Lipid metabolism": "#c7c7c7",
    "Terpenoid backbone biosynthesis": "#bcbd22",
    "Cofactor and vitamin metabolism": "#dbdb8d",
    "Polyamine biosynthesis": "#17becf",
    "ATP synthesis": "#9edae5",
    "Carbon fixation": "#393b79",
    "Sulfur metabolism": "#637939",
    "Methane metabolism": "#8c6d31",
    "Nitrogen metabolism": "#843c39",
    "Aromatics degradation": "#7b4173",
    "Pathogenicity": "#5254a3",
    "Drug resistance": "#6b6ecf",
    "Polyketide sugar unit biosynthesis": "#b5cf6b",
    "Biosynthesis of other antibiotics": "#9c9ede",
    "Other polysaccharide metabolism": "#cedb9c",
    "Other": "#d9d9d9",
}

CATEGORY_PATTERNS = [
    (category, re.compile(re.escape(category), re.IGNORECASE))
    for category in COLOR_MAP
    if category != "Other"
]

BACTERIA_PANELS = [
    ("All bacteria", lambda df: df["bacterium"] != "unclassified"),
    (
        "Carn. maltaromaticum",
        lambda df: df["bacterium"]
        == "g__Carnobacterium.s__Carnobacterium_maltaromaticum",
    ),
    (
        "Carn. divergens",
        lambda df: df["bacterium"]
        == "g__Carnobacterium.s__Carnobacterium_divergens",
    ),
    (
        "Vagoc. proximus",
        lambda df: df["bacterium"] == "g__Vagococcus.s__Vagococcus_proximus",
    ),
]


@dataclass(frozen=True)
class TemperatureConfig:
    input_csv: str
    output_png: str
    value_cols: tuple[str, ...]
    stage_columns: tuple[tuple[str, tuple[str, ...]], ...]
    interpolated_columns: tuple[tuple[str, str, str], ...] = ()
    acceptable_only_panels: tuple[str, ...] = ()
    acceptable_only_label: str | None = None


CONFIGS = {
    4: TemperatureConfig(
        input_csv="ko_all_reactions_4.csv",
        output_png="4_deg_bar_plots.png",
        value_cols=tuple(f"A{i:02d}" for i in range(1, 31)),
        stage_columns=(
            ("Acceptable \n Days 5-10", tuple(f"A{i:02d}" for i in range(1, 19))),
            ("Early spoiled \n Days 11-13", tuple(f"A{i:02d}" for i in range(19, 25))),
            ("Late spoiled \n Day 14", tuple(f"A{i:02d}" for i in range(25, 31))),
        ),
        interpolated_columns=(("A02", "A01", "A03"),),
    ),
    6: TemperatureConfig(
        input_csv="ko_all_reactions_6.csv",
        output_png="6_deg_bar_plots.png",
        value_cols=tuple(f"A{i:02d}" for i in range(5, 20)),
        stage_columns=(
            ("Acceptable \n Days 6-7", tuple(f"A{i:02d}" for i in range(5, 11))),
            ("Early spoiled \n Days 8-9", tuple(f"A{i:02d}" for i in range(11, 17))),
            ("Late spoiled \n Day 10", tuple(f"A{i:02d}" for i in range(17, 20))),
        ),
        # The 6 C Carn. maltaromaticum panel displays only the acceptable stage.
        acceptable_only_panels=("Carn. maltaromaticum",),
        acceptable_only_label="Acceptable \n Day 7",
    ),
}


def categorize_pathway(pathway_class: object) -> str:
    """Map a detailed pathway class to the categories shown in the legend."""
    if pd.isna(pathway_class):
        return "Other"
    pathway_class = str(pathway_class).split(" //")[0]
    for category, pattern in CATEGORY_PATTERNS:
        if pattern.search(pathway_class):
            return category
    return "Other"


def load_and_prepare(config: TemperatureConfig) -> pd.DataFrame:
    """Load reaction abundances and normalize every sample to TPM."""
    df = pd.read_csv(config.input_csv, index_col=0)
    for target, left, right in config.interpolated_columns:
        df[target] = (df[left] + df[right]) / 2
    df[list(config.value_cols)] = (
        df[list(config.value_cols)] / df[list(config.value_cols)].sum()
    ) * 1_000_000
    return df


def summarize_pathways(
    df: pd.DataFrame, config: TemperatureConfig
) -> pd.DataFrame:
    """Calculate each pathway category's proportion within each spoilage stage."""
    reactions = df[list(config.value_cols) + ["reaction"]].groupby("reaction").sum()
    reactions = pd.merge(
        reactions,
        df[["reaction", "class"]].drop_duplicates(),
        left_index=True,
        right_on="reaction",
    )
    reactions["Category"] = reactions["class"].apply(categorize_pathway)

    stage_values = {
        stage_label: reactions[list(columns)].mean(axis=1)
        for stage_label, columns in config.stage_columns
    }
    summary = pd.DataFrame(stage_values)
    summary["Category"] = reactions["Category"].values
    category_sums = summary.groupby("Category").sum()
    return category_sums.div(category_sums.sum(axis=0), axis=1)


def plot_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    config: TemperatureConfig,
    subtitle: str,
    panel_label: str,
    bar_width: float = 0.9,
    legend_fontsize: int = 11,
) -> None:
    """Plot one stacked pathway abundance panel."""
    proportions = summarize_pathways(df, config)
    stage_labels = [label for label, _ in config.stage_columns]
    if subtitle in config.acceptable_only_panels:
        stage_labels = stage_labels[:1]

    present_order = (
        proportions[stage_labels[0]]
        .loc[lambda values: values > 0]
        .sort_values(ascending=False)
        .index.tolist()
    )
    present_any = (
        proportions[stage_labels]
        .sum(axis=1)
        .loc[lambda values: values > 0]
        .index.tolist()
    )
    present_order += [category for category in present_any if category not in present_order]

    bottoms = {stage_label: 0.0 for stage_label in stage_labels}
    for category in present_order[::-1]:
        for stage_label in stage_labels:
            value = proportions.loc[category, stage_label]
            display_label = (
                config.acceptable_only_label
                if subtitle in config.acceptable_only_panels
                else stage_label
            )
            ax.bar(
                display_label,
                value,
                bottom=bottoms[stage_label],
                color=COLOR_MAP[category],
                edgecolor="white",
                width=bar_width,
            )
            bottoms[stage_label] += value

    ax.set_ylim(0, 1)
    ax.set_xlim(-0.5, len(stage_labels) - 0.5)
    ax.set_ylabel("Proportion of pathway transcripts abundance (TPM)")
    ax.set_title(subtitle)
    ax.text(
        0.02,
        1.1,
        panel_label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
    )
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLOR_MAP[category])
        for category in present_order
    ]
    ax.legend(
        handles,
        present_order,
        bbox_to_anchor=(1, 1),
        loc="upper left",
        fontsize=legend_fontsize,
    )


def create_pathway_figure(temperature: int) -> tuple[plt.Figure, Any]:
    """Create and save the final a)-d) pathway plot for 4 C or 6 C."""
    try:
        config = CONFIGS[temperature]
    except KeyError as exc:
        raise ValueError(f"Unsupported temperature: {temperature}. Use 4 or 6.") from exc

    df = load_and_prepare(config)
    fig, axes = plt.subplots(2, 2, figsize=(16, 15.5))
    fig.suptitle(
        f"Functional activity distribution {temperature} C: "
        "Acceptable vs Early spoiled vs Late spoiled",
        fontsize=14,
    )

    for ax, panel_label, (subtitle, filter_expr) in zip(
        axes.flat, "abcd", BACTERIA_PANELS
    ):
        plot_panel(
            ax,
            df.loc[filter_expr(df)],
            config,
            subtitle,
            panel_label,
        )

    fig.tight_layout(rect=[0, 0, 0.92, 0.95])
    fig.savefig(config.output_png)
    plt.show()
    return fig, axes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--temperature",
        type=int,
        choices=sorted(CONFIGS),
        required=True,
        help="Storage temperature in degrees Celsius.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    create_pathway_figure(args.temperature)
