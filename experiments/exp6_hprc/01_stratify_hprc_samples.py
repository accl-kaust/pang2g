#!/usr/bin/env python3
"""
01_stratify_hprc_samples.py — Assign HPRC samples to AFR/EAS populations.

Classification is taken STRICTLY from the authoritative 1000 Genomes
3202-sample population table (downloaded by download_hprc_assemblies.sh to
$G2G_DATA_ROOT/hprc/metadata/1kgp_pops.panel) — never from a hardcoded guess.
An earlier hardcoded map mislabeled several American samples (HG00733=PUR,
HG01123/HG01175=CLM, HG02148=PEL) as EAS, which would have built the "EAS"
pangenome from the wrong genomes and invalidated the AFR-vs-EAS comparison.

For each population we select an EQUAL number of samples (balanced, capped by
availability and --n-per-pop), keeping only samples whose BOTH haplotype
assemblies are present on disk. The selection matches what
download_hprc_assemblies.sh already wrote; running this is idempotent and acts
as a verification/regeneration step on the compute node.

Usage:
    python 01_stratify_hprc_samples.py \\
        --panel        ~/g2g_data/hprc/metadata/1kgp_pops.panel \\
        --index        ~/g2g_data/hprc/metadata/hprc_v2_index.tsv \\
        --assembly-dir ~/g2g_data/hprc/assemblies/ \\
        --out-dir      ~/g2g_data/hprc/metadata/
"""

import argparse
import csv
import os
import sys

POPS = ("AFR", "EAS")


def load_superpop(panel_file: str) -> dict[str, str]:
    """sample_id -> super-population, parsed from the 1KGP table by header name.

    Handles both the 3202 ped format (space-separated:
    FamilyID SampleID FatherID MotherID Sex Population Superpopulation) and the
    2504 panel (tab-separated: sample pop super_pop gender).
    """
    superpop: dict[str, str] = {}
    if not os.path.exists(panel_file):
        print(f"[stratify] ERROR: 1KGP panel not found: {panel_file}", file=sys.stderr)
        print("[stratify] Run download_hprc_assemblies.sh first.", file=sys.stderr)
        sys.exit(2)

    with open(panel_file) as f:
        first = f.readline()
        sep = "\t" if "\t" in first else None  # None => any whitespace
        hdr = [h.lower() for h in first.split(sep)]

        def col(*names):
            for i, h in enumerate(hdr):
                if h in names:
                    return i
            return None

        s_i = col("sample", "sampleid", "sample_id")
        sp_i = col("super_pop", "superpopulation", "superpop", "super_population")
        if s_i is None or sp_i is None:
            print(f"[stratify] ERROR: cannot find sample/superpop columns in "
                  f"panel header: {hdr}", file=sys.stderr)
            sys.exit(2)
        for line in f:
            p = line.split(sep)
            if len(p) > max(s_i, sp_i):
                superpop[p[s_i].strip()] = p[sp_i].strip()

    # GIAB Han Chinese son (NA24631) — not present in 1KGP tables.
    superpop.setdefault("HG005", "EAS")
    print(f"[stratify] super-pop entries loaded: {len(superpop)}", file=sys.stderr)
    return superpop


def load_index_samples(index_file: str) -> list[str]:
    """Sample IDs present in the HPRC assembly index."""
    if not os.path.exists(index_file):
        print(f"[stratify] ERROR: index not found: {index_file}", file=sys.stderr)
        sys.exit(2)
    samples = []
    with open(index_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        scol = next((c for c in (reader.fieldnames or [])
                     if c.lower() in ("sample", "sample_id", "sampleid")),
                    (reader.fieldnames or [None])[0])
        for row in reader:
            s = (row.get(scol) or "").strip()
            if s:
                samples.append(s)
    return samples


def both_haplotypes_present(assembly_dir: str, pop: str, sample: str) -> bool:
    h1 = os.path.join(assembly_dir, pop, sample, f"{sample}_hap1.fa.gz")
    h2 = os.path.join(assembly_dir, pop, sample, f"{sample}_hap2.fa.gz")
    return os.path.exists(h1) and os.path.exists(h2)


def stratify(panel_file: str, index_file: str, assembly_dir: str,
             out_dir: str, n_per_pop: int) -> dict[str, list[str]]:
    os.makedirs(out_dir, exist_ok=True)
    superpop = load_superpop(panel_file)
    index_samples = load_index_samples(index_file)

    avail = {p: [] for p in POPS}
    for s in index_samples:
        sp = superpop.get(s)
        if sp in avail:
            if assembly_dir and not both_haplotypes_present(assembly_dir, sp, s):
                print(f"[stratify] {s} ({sp}): missing haplotype(s), skipping",
                      file=sys.stderr)
                continue
            avail[sp].append(s)

    n = min(min(len(avail[p]) for p in POPS), n_per_pop)
    print(f"[stratify] available " +
          " ".join(f"{p}={len(avail[p])}" for p in POPS) +
          f" -> balanced n={n} per pop", file=sys.stderr)
    if n == 0:
        print("[stratify] ERROR: could not select a balanced AFR/EAS set",
              file=sys.stderr)
        sys.exit(3)

    result = {p: sorted(avail[p])[:n] for p in POPS}

    for pop, samples in result.items():
        out_file = os.path.join(out_dir, f"{pop.lower()}_samples.txt")
        with open(out_file, "w") as f:
            f.write("\n".join(samples) + "\n")
        print(f"[stratify] {pop}: {len(samples)} samples -> {out_file}")

    combined = os.path.join(out_dir, "all_samples_with_pop.tsv")
    with open(combined, "w") as f:
        f.write("sample_id\tpopulation\n")
        for pop, samples in result.items():
            for s in samples:
                f.write(f"{s}\t{pop}\n")
    print(f"[stratify] Combined table -> {combined}")
    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    data_root = os.environ.get("G2G_DATA_ROOT", os.path.expanduser("~/g2g_data"))
    md = os.path.join(data_root, "hprc/metadata")
    parser.add_argument("--panel", default=os.path.join(md, "1kgp_pops.panel"))
    parser.add_argument("--index", default=os.path.join(md, "hprc_v2_index.tsv"))
    parser.add_argument("--assembly-dir",
                        default=os.path.join(data_root, "hprc/assemblies/"))
    parser.add_argument("--out-dir", default=md + "/")
    parser.add_argument("--n-per-pop", type=int,
                        default=int(os.environ.get("N_PER_POP", "20")))
    # Back-compat: accept the old --metadata flag as an alias for --index.
    parser.add_argument("--metadata", dest="index_compat", default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--n-afr", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--n-eas", type=int, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    index = args.index_compat or args.index
    n_per_pop = args.n_per_pop
    if args.n_afr is not None or args.n_eas is not None:
        n_per_pop = min(x for x in (args.n_afr, args.n_eas, n_per_pop)
                        if x is not None)

    result = stratify(
        panel_file=args.panel,
        index_file=index,
        assembly_dir=args.assembly_dir,
        out_dir=args.out_dir,
        n_per_pop=n_per_pop,
    )

    print("\n[stratify] Final counts:")
    for pop, s in result.items():
        print(f"  {pop}: {len(s)} samples  {s}")


if __name__ == "__main__":
    main()
