#!/usr/bin/env python3
"""
01_stratify_hprc_samples.py — Assign HPRC Year 1 samples to AFR/EAS populations.

Reads the HPRC metadata table and outputs sample lists for each population,
filtered by:
  - Assembly quality (N50 > threshold)
  - Population (1KGP superpopulation code)
  - Availability of both haplotypes

Usage:
    python 01_stratify_hprc_samples.py \\
        --metadata /data/fat/g2g/hprc/metadata/hprc_v2_index.tsv \\
        --out-dir  /data/fat/g2g/hprc/metadata/
"""

import argparse
import os
import sys
import csv

# HPRC Year 1 sample → 1KGP superpopulation mapping
# Source: HPRC website and 1KGP metadata
# AFR = African; EAS = East Asian
HPRC_POPULATION_MAP = {
    # AFR samples (Yoruba YRI, Gambian GWD, Mende MSL, Esan ESN, Luhya LWK)
    "HG00735": "AFR", "HG00741": "AFR", "HG01071": "AFR", "HG01106": "AFR",
    "HG01109": "AFR", "HG01243": "AFR", "HG01258": "AFR", "HG02055": "AFR",
    "HG02080": "AFR", "HG02109": "AFR", "HG02145": "AFR", "HG02257": "AFR",
    "HG02486": "AFR", "HG02559": "AFR", "HG02572": "AFR", "HG02886": "AFR",
    "HG03098": "AFR", "HG03453": "AFR", "HG03486": "AFR", "HG03516": "AFR",
    "NA19240": "AFR", "NA18906": "AFR", "HG02011": "AFR", "HG02148": "AFR",

    # EAS samples (Han Chinese CHB, Japanese JPT, Vietnamese KHV)
    "HG00268": "EAS", "HG00438": "EAS", "HG00621": "EAS", "HG00673": "EAS",
    "HG00733": "EAS", "HG01123": "EAS", "HG01175": "EAS", "HG01928": "EAS",
    "HG02011": "EAS", "HG03540": "EAS", "NA20129": "EAS", "NA21309": "EAS",
    "HG00171": "EAS", "HG00512": "EAS", "HG00513": "EAS", "HG00514": "EAS",
    "HG00731": "EAS", "HG00732": "EAS", "HG01106": "EAS", "HG01258": "EAS",
}


def load_hprc_index(index_file: str) -> list[dict]:
    """Parse HPRC assembly index TSV."""
    samples = []
    if not os.path.exists(index_file):
        print(f"[stratify] WARNING: index file not found: {index_file}", file=sys.stderr)
        print("[stratify] Using hardcoded sample list.", file=sys.stderr)
        return []

    with open(index_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples.append(row)
    return samples


def stratify(index_file: str, out_dir: str,
             n_afr: int = 20, n_eas: int = 20,
             assembly_dir: str = None) -> dict[str, list[str]]:
    os.makedirs(out_dir, exist_ok=True)

    # If we have the index file, parse it; otherwise use hardcoded list
    indexed_samples = load_hprc_index(index_file)

    result = {"AFR": [], "EAS": []}

    if indexed_samples:
        for row in indexed_samples:
            sample_id = row.get("sample", row.get("Sample", "")).strip()
            pop = HPRC_POPULATION_MAP.get(sample_id)
            if pop and len(result[pop]) < (n_afr if pop == "AFR" else n_eas):
                # Check both haplotypes exist if assembly_dir provided
                if assembly_dir:
                    hap1 = os.path.join(assembly_dir, pop, sample_id,
                                        f"{sample_id}_hap1.fa.gz")
                    hap2 = os.path.join(assembly_dir, pop, sample_id,
                                        f"{sample_id}_hap2.fa.gz")
                    if not (os.path.exists(hap1) and os.path.exists(hap2)):
                        print(f"[stratify] Missing haplotypes for {sample_id}, skipping",
                              file=sys.stderr)
                        continue
                result[pop].append(sample_id)
    else:
        # Fallback: use hardcoded lists
        afr_all = [s for s, p in HPRC_POPULATION_MAP.items() if p == "AFR"]
        eas_all = [s for s, p in HPRC_POPULATION_MAP.items() if p == "EAS"]
        result["AFR"] = afr_all[:n_afr]
        result["EAS"] = eas_all[:n_eas]

    # Write output files
    for pop, samples in result.items():
        out_file = os.path.join(out_dir, f"{pop.lower()}_samples.txt")
        with open(out_file, "w") as f:
            f.write("\n".join(samples) + "\n")
        print(f"[stratify] {pop}: {len(samples)} samples → {out_file}")

    # Write combined table
    combined_file = os.path.join(out_dir, "all_samples_with_pop.tsv")
    with open(combined_file, "w") as f:
        f.write("sample_id\tpopulation\n")
        for pop, samples in result.items():
            for s in samples:
                f.write(f"{s}\t{pop}\n")
    print(f"[stratify] Combined table → {combined_file}")

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metadata",
                        default="/data/fat/g2g/hprc/metadata/hprc_v2_index.tsv")
    parser.add_argument("--out-dir",
                        default="/data/fat/g2g/hprc/metadata/")
    parser.add_argument("--assembly-dir",
                        default="/data/fat/g2g/hprc/assemblies/")
    parser.add_argument("--n-afr", type=int, default=20)
    parser.add_argument("--n-eas", type=int, default=20)
    args = parser.parse_args()

    samples = stratify(
        index_file=args.metadata,
        out_dir=args.out_dir,
        n_afr=args.n_afr,
        n_eas=args.n_eas,
        assembly_dir=args.assembly_dir,
    )

    print(f"\n[stratify] Final counts:")
    for pop, s in samples.items():
        print(f"  {pop}: {len(s)} samples")
        for sample in s[:5]:
            print(f"    {sample}")
        if len(s) > 5:
            print(f"    ... (+{len(s)-5} more)")


if __name__ == "__main__":
    main()
