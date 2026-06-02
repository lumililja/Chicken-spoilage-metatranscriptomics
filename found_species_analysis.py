"""Create MetaPhlAn-style profiles from Bracken output and check ChocoPhlAn."""

from pathlib import Path

import pandas as pd


BRACKEN_DIR = Path("bracken_outputs_6")
OUTPUT_DIR = Path("bracken_to_metaphlan")
CHOCOPHLAN_DIR = Path("../HUMANN/HUMANN_NUC")
MINIMUM_ABUNDANCE = 0.05
PROFILE_HEADER = "#mpa_vJun23_CHOCOPhlAnSGB_202403"

SPECIES_RENAMES = {
    "Lactococcus_carnosus": "Lactococcus_piscium",
    "Latilactobacillus_sakei": "Lactobacillus_sakei",
}


def create_profiles() -> set[str]:
    """Write filtered MetaPhlAn-style files and return all identified species."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    all_species = set()

    for input_file in sorted(BRACKEN_DIR.glob("*.bracken")):
        df = pd.read_csv(input_file, sep="\t")
        df = df.loc[df["fraction_total_reads"] >= MINIMUM_ABUNDANCE].copy()
        df["name"] = df["name"].str.replace(" ", "_", regex=False)
        df["name"] = df["name"].replace(SPECIES_RENAMES)
        all_species.update(f"s__{species}" for species in df["name"])
        df["name"] = "|s__" + df["name"]

        output_file = OUTPUT_DIR / f"{input_file.stem}_metaphlan.txt"
        with output_file.open("w") as profile:
            profile.write(f"{PROFILE_HEADER}\n")
            df[["name", "taxonomy_id", "fraction_total_reads", "new_est_reads"]].to_csv(
                profile, sep="\t", header=False, index=False
            )

    return all_species


def report_chocophlan_species(species: set[str]) -> None:
    """Print species found and missing from the ChocoPhlAn nucleotide database."""
    if not CHOCOPHLAN_DIR.is_dir():
        raise FileNotFoundError(f"ChocoPhlAn database directory not found: {CHOCOPHLAN_DIR}")

    genome_names = [genome.name for genome in CHOCOPHLAN_DIR.iterdir()]
    found = sorted(s for s in species if any(s in genome for genome in genome_names))
    missing = sorted(species - set(found))

    print("\nSpecies found in ChocoPhlAn:")
    print("\n".join(found) or "None")
    print("\nSpecies missing from ChocoPhlAn:")
    print("\n".join(missing) or "None")


if __name__ == "__main__":
    identified_species = create_profiles()
    print(f"Wrote {len(list(OUTPUT_DIR.glob('*_metaphlan.txt')))} profiles.")
    report_chocophlan_species(identified_species)
